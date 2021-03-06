
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// integrators/sppm.cpp*
#include "integrators/vsppm.h"
#include "scene.h"
#include "imageio.h"
#include "paramset.h"
#include "progressreporter.h"
#include "sampling.h"
#include "samplers/halton.h"

namespace pbrt {

STAT_RATIO(
        "Stochastic Progressive Photon Mapping/Visible points checked per photon "
                "intersection",
        visiblePointsChecked, totalPhotonSurfaceInteractions);
STAT_COUNTER("Stochastic Progressive Photon Mapping/Photon paths followed",
             photonPaths);
STAT_COUNTER("Stochastic Progressive Photon Mapping/Total photon medium interactions",
             totalPhotonMediumInteractions);
STAT_COUNTER("Stochastic Progressive Photon Mapping/Total visible points in participating media",
             totalVisibleMediumPoints);
STAT_COUNTER("Stochastic Progressive Photon Mapping/Total visible points on surfaces",
             totalVisibleSurfacePoints);
STAT_INT_DISTRIBUTION(
        "Stochastic Progressive Photon Mapping/Grid cells per visible point",
        gridCellsPerVisiblePoint);
STAT_MEMORY_COUNTER("Memory/SPPM Pixels", pixelMemoryBytes);
STAT_FLOAT_DISTRIBUTION("Memory/SPPM BSDF and Grid Memory", memoryArenaMB);

struct VisiblePoint {
    // VisiblePoint Public Methods
    VisiblePoint() : bsdf(nullptr), phase(nullptr) {}

    Point3f p;
    Vector3f wo;
    const BSDF *bsdf;
    const PhaseFunction *phase;
    Spectrum beta;

    bool IsSurfacePoint() const {
        if (bsdf) CHECK(!phase);
        return bsdf != nullptr;
    }

    bool IsMediumPoint() const {
        if (phase) CHECK(!bsdf);
        return phase != nullptr;
    }
};

// SPPM Local Definitions
struct SPPMPixel {
    // SPPMPixel Public Methods
    SPPMPixel() : M(0) {}

    // SPPMPixel Public Data
    Float radius = 0;
    Spectrum Ld;
    VisiblePoint vp;
    AtomicFloat Phi[Spectrum::nSamples];
    std::atomic<int> M;
    Float N = 0;
    Spectrum tau;
};

struct SPPMPixelListNode {
    SPPMPixel *pixel;
    SPPMPixelListNode *next;
};

static bool ToGrid(const Point3f &p, const Bounds3f &bounds,
                   const int gridRes[3], Point3i *pi) {
    bool inBounds = true;
    Vector3f pg = bounds.Offset(p);
    for (int i = 0; i < 3; ++i) {
        (*pi)[i] = (int) (gridRes[i] * pg[i]);
        inBounds &= ((*pi)[i] >= 0 && (*pi)[i] < gridRes[i]);
        (*pi)[i] = Clamp((*pi)[i], 0, gridRes[i] - 1);
    }
    return inBounds;
}

inline unsigned int hash(const Point3i &p, int hashSize) {
    return (unsigned int) ((p.x * 73856093) ^ (p.y * 19349663) ^
                           (p.z * 83492791)) %
           hashSize;
}

class AwesomeSampler : public Sampler {
public:
    AwesomeSampler(int64_t samplesPerPixel,
                   Sampler &sampler,
                   size_t const SamplerLimit,
                   uint64_t rngSeed) : Sampler(samplesPerPixel),
                                       rng(rngSeed),
                                       sampler(sampler),
                                       SamplerLimit(SamplerLimit),
                                       sampleCount(0) {}

    Float Get1D() override {
        sampleCount++;
        if (sampleCount <= SamplerLimit) {
            return sampler.Get1D();
        }
        return rng.UniformFloat();
    }

    Point2f Get2D() override {
        sampleCount += 2;
        if (sampleCount <= SamplerLimit) {
            return sampler.Get2D();
        }
        return Point2f(rng.UniformFloat(), rng.UniformFloat());
    }

    std::unique_ptr<Sampler> Clone(int seed) override {
        return nullptr;
    }

private:
    RNG rng;
    Sampler &sampler;
    const size_t SamplerLimit;
    size_t sampleCount;
};

class AwesomeHaltonSampler : public Sampler {
public:
    AwesomeHaltonSampler(uint64_t haltonIndex)
            : Sampler(0), HaltonIndex(haltonIndex), haltonDim(0), rng(HaltonIndex) {}

    Float Get1D() override {
        if (haltonDim + 1 <= 1000) {
            return RadicalInverse(haltonDim++, HaltonIndex);
        }
        return rng.UniformFloat();
    }

    Point2f Get2D() override {
        return Point2f(Get1D(), Get1D());
    }

    std::unique_ptr<Sampler> Clone(int seed) override {
        return nullptr;
    }

private:
    const uint64_t HaltonIndex;
    uint32_t haltonDim;
    RNG rng;
};

// SPPM Method Definitions
void VolSPPMIntegrator::Render(const Scene &scene) {
    ProfilePhase p(Prof::IntegratorRender);
    // Initialize _pixelBounds_ and _pixels_ array for SPPM
    Bounds2i pixelBounds = camera->film->croppedPixelBounds;
    int nPixels = pixelBounds.Area();
    std::unique_ptr<SPPMPixel[]> pixels(new SPPMPixel[nPixels]);
    for (int i = 0; i < nPixels; ++i) pixels[i].radius = initialSearchRadius;
    const Float invSqrtSPP = 1.f / std::sqrt(nIterations);
    pixelMemoryBytes = nPixels * sizeof(SPPMPixel);
    // Compute _lightDistr_ for sampling lights proportional to power
    std::unique_ptr<Distribution1D> lightDistr =
            ComputeLightPowerDistribution(scene);

    // Perform _nIterations_ of VolSPPM integration
    HaltonSampler sampler(nIterations, pixelBounds);

    // Compute number of tiles to use for VolSPPM camera pass
    Vector2i pixelExtent = pixelBounds.Diagonal();
    const int tileSize = 16;
    Point2i nTiles((pixelExtent.x + tileSize - 1) / tileSize,
                   (pixelExtent.y + tileSize - 1) / tileSize);
    ProgressReporter progress(2 * nIterations, "Rendering");

    // Counter to extend HaltonSampler values beyond 1000
    uint64_t globalNumPixels = 0;

    for (int iter = 0; iter < nIterations; ++iter) {
        // Counter to extend HaltonSampler values beyond 1000
        uint32_t iterNumPixels = 0;

        // Generate VolSPPM visible points
        std::vector<MemoryArena> perThreadArenas(MaxThreadIndex());
        {
            ProfilePhase _(Prof::SPPMCameraPass);
            ParallelFor2D([&](Point2i tile) {
                MemoryArena &arena = perThreadArenas[ThreadIndex];
                // Follow camera paths for _tile_ in image for VolSPPM
                const int TileIndex = tile.y * nTiles.x + tile.x;
                std::unique_ptr<Sampler> tileSampler = sampler.Clone(TileIndex);
                // Compute _tileBounds_ for VolSPPM tile
                int x0 = pixelBounds.pMin.x + tile.x * tileSize;
                int x1 = std::min(x0 + tileSize, pixelBounds.pMax.x);
                int y0 = pixelBounds.pMin.y + tile.y * tileSize;
                int y1 = std::min(y0 + tileSize, pixelBounds.pMax.y);
                Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));

                for (Point2i pPixel : tileBounds) {
                    const uint32_t PixelIndex = iterNumPixels++;
                    const uint64_t GoodPixelIndex = globalNumPixels + PixelIndex;
                    // Prepare _tileSampler_ for _pPixel_
                    tileSampler->StartPixel(pPixel);
                    tileSampler->SetSampleNumber(iter);
                    AwesomeSampler localSampler(0, *tileSampler, 1000, GoodPixelIndex);

                    // Generate camera ray for pixel for SPPM
                    CameraSample cameraSample =
                            localSampler.GetCameraSample(pPixel);
                    RayDifferential ray;
                    Spectrum beta =
                            camera->GenerateRayDifferential(cameraSample, &ray);
                    ray.ScaleDifferentials(invSqrtSPP);

                    // Follow camera ray path until a visible point is created

                    // Get _SPPMPixel_ for _pPixel_
                    Point2i pPixelO = Point2i(pPixel - pixelBounds.pMin);
                    int pixelOffset =
                            pPixelO.x +
                            pPixelO.y * (pixelBounds.pMax.x - pixelBounds.pMin.x);
                    SPPMPixel &pixel = pixels[pixelOffset];
                    bool specularBounce = false;
                    for (int depth = 0; depth < maxDepth; ++depth) {
                        SurfaceInteraction isect;
                        if (!scene.Intersect(ray, &isect)) {
                            // Accumulate light contributions for ray with no
                            // intersection
                            for (const auto &light : scene.lights)
                                pixel.Ld += beta * light->Le(ray);
                            break;
                        }

                        MediumInteraction mi;
                        if (ray.medium) {
                            if (renderMedia)
                                beta *= ray.medium->Sample(ray, localSampler, arena, &mi);
                            else
                                beta *= ray.medium->Tr(ray, localSampler);
                        }

                        if (beta.IsBlack()) break;

                        // Handle an interaction with a medium or a surface
                        if (renderMedia && mi.IsValid()) {
                            ++totalPhotonMediumInteractions;
                            // Process SPPM camera ray medium interaction
                            Vector3f wo = -ray.d;
                            pixel.Ld += beta * UniformSampleOneLight(mi, scene, arena,
                                                                     localSampler, true);

                            // Create visible point and end camera path
                            pixel.vp.p = mi.p;
                            pixel.vp.wo = wo;
                            pixel.vp.bsdf = nullptr;
                            pixel.vp.phase = mi.phase;
                            pixel.vp.beta = beta;
                            ++totalVisibleMediumPoints;
                            break;

                        } else if (renderSurfaces) {
                            ++totalPhotonSurfaceInteractions;
                            // Process SPPM camera ray surface intersection

                            // Compute BSDF at SPPM camera ray intersection
                            isect.ComputeScatteringFunctions(ray, arena, true);
                            if (!isect.bsdf) {
                                ray = isect.SpawnRay(ray.d);
                                --depth;
                                continue;
                            }
                            const BSDF &bsdf = *isect.bsdf;

                            // Accumulate direct illumination at SPPM camera ray
                            // intersection
                            Vector3f wo = -ray.d;
                            if (depth == 0 || specularBounce)
                                pixel.Ld += beta * isect.Le(wo);
                            pixel.Ld +=
                                    beta * UniformSampleOneLight(isect, scene, arena,
                                                                 localSampler, true);

                            // Possibly create visible point and end camera path
                            bool isDiffuse = bsdf.NumComponents(BxDFType(
                                    BSDF_DIFFUSE | BSDF_REFLECTION |
                                    BSDF_TRANSMISSION)) > 0;
                            bool isGlossy = bsdf.NumComponents(BxDFType(
                                    BSDF_GLOSSY | BSDF_REFLECTION |
                                    BSDF_TRANSMISSION)) > 0;
                            if (isDiffuse || (isGlossy && depth == maxDepth - 1)) {
                                pixel.vp.p = isect.p;
                                pixel.vp.wo = wo;
                                pixel.vp.bsdf = isect.bsdf;
                                pixel.vp.phase = nullptr;
                                pixel.vp.beta = beta;
                                ++totalVisibleSurfacePoints;
                                break;
                            }

                            // Spawn ray from SPPM camera path vertex
                            if (depth < maxDepth - 1) {
                                Float pdf;
                                Vector3f wi;
                                BxDFType type;
                                Spectrum f =
                                        bsdf.Sample_f(wo, &wi, localSampler.Get2D(),
                                                      &pdf, BSDF_ALL, &type);
                                if (pdf == 0. || f.IsBlack()) break;
                                specularBounce = (type & BSDF_SPECULAR) != 0;
                                beta *= f * AbsDot(wi, isect.shading.n) / pdf;
                                if (beta.y() < 0.25) {
                                    Float continueProb =
                                            std::min((Float) 1, beta.y());
                                    if (localSampler.Get1D() > continueProb) break;
                                    beta /= continueProb;
                                }
                                ray = (RayDifferential) isect.SpawnRay(wi);
                            }
                        }
                    }
                }
            }, nTiles);
        }
        progress.Update();

        // Create grid of all SPPM visible points
        int gridRes[3];
        Bounds3f gridBounds;
        // Allocate grid for SPPM visible points
        const int hashSize = nPixels;
        std::vector<std::atomic<SPPMPixelListNode *>> grid(hashSize);
        {
            ProfilePhase _(Prof::SPPMGridConstruction);

            // Compute grid bounds for SPPM visible points
            Float maxRadius = 0.;
            for (int i = 0; i < nPixels; ++i) {
                const SPPMPixel &pixel = pixels[i];
                if (pixel.vp.beta.IsBlack()) continue;
                Bounds3f vpBound = Expand(Bounds3f(pixel.vp.p), pixel.radius);
                gridBounds = Union(gridBounds, vpBound);
                maxRadius = std::max(maxRadius, pixel.radius);
            }

            // Compute resolution of SPPM grid in each dimension
            const Vector3f diag = gridBounds.Diagonal();
            const Float maxDiag = MaxComponent(diag);
            const int baseGridRes = (int) (maxDiag / maxRadius);
            CHECK_GT(baseGridRes, 0);
            for (int i = 0; i < 3; ++i)
                gridRes[i] = std::max((int) (baseGridRes * diag[i] / maxDiag), 1);

            // Add visible points to SPPM grid
            ParallelFor([&](int pixelIndex) {
                MemoryArena &arena = perThreadArenas[ThreadIndex];
                SPPMPixel &pixel = pixels[pixelIndex];
                if (!pixel.vp.beta.IsBlack()) {
                    // Add pixel's visible point to applicable grid cells
                    Float radius = pixel.radius;
                    Point3i pMin, pMax;
                    ToGrid(pixel.vp.p - Vector3f(radius, radius, radius),
                           gridBounds, gridRes, &pMin);
                    ToGrid(pixel.vp.p + Vector3f(radius, radius, radius),
                           gridBounds, gridRes, &pMax);
                    for (int z = pMin.z; z <= pMax.z; ++z)
                        for (int y = pMin.y; y <= pMax.y; ++y)
                            for (int x = pMin.x; x <= pMax.x; ++x) {
                                // Add visible point to grid cell $(x, y, z)$
                                int h = hash(Point3i(x, y, z), hashSize);
                                SPPMPixelListNode *node =
                                        arena.Alloc<SPPMPixelListNode>();
                                node->pixel = &pixel;

                                // Atomically add _node_ to the start of
                                // _grid[h]_'s linked list
                                node->next = grid[h];
                                while (!grid[h].compare_exchange_weak(node->next, node));
                            }
                    ReportValue(gridCellsPerVisiblePoint,
                                (1 + pMax.x - pMin.x) * (1 + pMax.y - pMin.y) *
                                (1 + pMax.z - pMin.z));
                }
            }, nPixels, 4096);
        }

        // Trace photons and accumulate contributions
        {
            ProfilePhase _(Prof::SPPMPhotonPass);
            std::vector<MemoryArena> photonShootArenas(MaxThreadIndex());
            ParallelFor([&](int photonIndex) {
                MemoryArena &arena = photonShootArenas[ThreadIndex];
                // Follow photon path for _photonIndex_
                const uint64_t HaltonIndex =
                        (uint64_t) iter * (uint64_t) photonsPerIteration +
                        photonIndex;
                AwesomeHaltonSampler localSampler(HaltonIndex);

                // Choose light to shoot photon from
                Float lightPdf;
                Float lightSample = localSampler.Get1D();
                int lightNum =
                        lightDistr->SampleDiscrete(lightSample, &lightPdf);
                const std::shared_ptr<Light> &light = scene.lights[lightNum];

                // Compute sample values for photon ray leaving light source
                Point2f uLight0 = localSampler.Get2D();
                Point2f uLight1 = localSampler.Get2D();
                Float uLightTime = localSampler.Get1D();

                // Generate _photonRay_ from light source and initialize _beta_
                RayDifferential photonRay;
                Normal3f nLight;
                Float pdfPos, pdfDir;
                Spectrum Le =
                        light->Sample_Le(uLight0, uLight1, uLightTime, &photonRay,
                                         &nLight, &pdfPos, &pdfDir);
                if (pdfPos == 0 || pdfDir == 0 || Le.IsBlack()) return;
                Spectrum beta = (AbsDot(nLight, photonRay.d) * Le) /
                                (lightPdf * pdfPos * pdfDir);
                if (beta.IsBlack()) return;

                // Follow photon path through scene and record intersections
                SurfaceInteraction isect;
                for (int depth = 0; depth < maxDepth; ++depth) {
                    if (!scene.Intersect(photonRay, &isect)) break;

                    MediumInteraction mi;
                    if (photonRay.medium) beta *= photonRay.medium->Sample(photonRay, localSampler, arena, &mi);
                    if (beta.IsBlack()) break;

                    // Handle an interaction with a medium or a surface
                    Spectrum bnew;
                    Vector3f newWi;
                    if (mi.IsValid()) {
                        ++totalPhotonMediumInteractions;
                        Point3i photonGridIndex;
                        if (ToGrid(mi.p, gridBounds, gridRes,
                                   &photonGridIndex)) {
                            int h = hash(photonGridIndex, hashSize);
                            // Add photon contribution to visible points in
                            // _grid[h]_
                            for (SPPMPixelListNode *node =
                                    grid[h].load(std::memory_order_relaxed);
                                 node != nullptr; node = node->next) {
                                ++visiblePointsChecked;
                                SPPMPixel &pixel = *node->pixel;
                                if (pixel.vp.IsMediumPoint()) {
                                    Float radius = pixel.radius;
                                    if (DistanceSquared(pixel.vp.p, mi.p) >
                                        radius * radius)
                                        continue;
                                    // Update _pixel_ $\Phi$ and $M$ for nearby
                                    // photon
                                    Vector3f wi = -photonRay.d;
                                    Spectrum Phi = beta * pixel.vp.phase->p(pixel.vp.wo, wi);
                                    for (int i = 0; i < Spectrum::nSamples; ++i)
                                        pixel.Phi[i].Add(Phi[i]);
                                    ++pixel.M;
                                }
                            }

                            // Sample new photon ray direction from volume
                            Vector3f wo = -photonRay.d;
                            mi.phase->Sample_p(wo, &newWi, localSampler.Get2D());
                            photonRay = mi.SpawnRay(newWi);
                        }

                    } else {
                        ++totalPhotonSurfaceInteractions;
                        if (depth > 0) {
                            // Add photon contribution to nearby visible points
                            Point3i photonGridIndex;
                            if (ToGrid(isect.p, gridBounds, gridRes,
                                       &photonGridIndex)) {
                                int h = hash(photonGridIndex, hashSize);
                                // Add photon contribution to visible points in
                                // _grid[h]_
                                for (SPPMPixelListNode *node =
                                        grid[h].load(std::memory_order_relaxed);
                                     node != nullptr; node = node->next) {
                                    ++visiblePointsChecked;
                                    SPPMPixel &pixel = *node->pixel;
                                    if (pixel.vp.IsSurfacePoint()) {
                                        Float radius = pixel.radius;
                                        if (DistanceSquared(pixel.vp.p, isect.p) >
                                            radius * radius)
                                            continue;
                                        // Update _pixel_ $\Phi$ and $M$ for nearby
                                        // photon
                                        Vector3f wi = -photonRay.d;
                                        Spectrum Phi =
                                                beta * pixel.vp.bsdf->f(pixel.vp.wo, wi);
                                        for (int i = 0; i < Spectrum::nSamples; ++i)
                                            pixel.Phi[i].Add(Phi[i]);
                                        ++pixel.M;
                                    }
                                }
                            }
                        }
                        // Sample new photon ray direction from surface

                        // Compute BSDF at photon intersection point
                        isect.ComputeScatteringFunctions(photonRay, arena, true,
                                                         TransportMode::Importance);
                        if (!isect.bsdf) {
                            --depth;
                            photonRay = isect.SpawnRay(photonRay.d);
                            continue;
                        }
                        const BSDF &photonBSDF = *isect.bsdf;

                        // Sample BSDF _fr_ and direction _wi_ for reflected photon
                        Vector3f wo = -photonRay.d;
                        Float pdf;
                        BxDFType flags;

                        // Generate _bsdfSample_ for outgoing photon sample
                        Point2f bsdfSample = localSampler.Get2D();
                        Spectrum fr = photonBSDF.Sample_f(wo, &newWi, bsdfSample, &pdf,
                                                          BSDF_ALL, &flags);
                        if (fr.IsBlack() || pdf == 0.f) break;
                        bnew = beta * fr * AbsDot(newWi, isect.shading.n) / pdf;

                        photonRay = (RayDifferential) isect.SpawnRay(newWi);
                    }

                    // Possibly terminate photon path with Russian roulette
                    Float q = std::max((Float) 0, 1 - bnew.y() / beta.y());
                    if (localSampler.Get1D() < q) break;
                    beta = bnew / (1 - q);
                }
                arena.Reset();
            }, photonsPerIteration, 8192);
            progress.Update();
            photonPaths += photonsPerIteration;
        }

        // Update pixel values from this pass's photons
        {
            ProfilePhase _(Prof::SPPMStatsUpdate);
            ParallelFor([&](int i) {
                SPPMPixel &p = pixels[i];
                if (((p.vp.IsSurfacePoint() && renderSurfaces) ||
                    (p.vp.IsMediumPoint() && renderMedia)) && p.M > 0) {
                    // Update pixel photon count, search radius, and $\tau$ from
                    // photons
                    Float gamma = (Float) 2 / (Float) 3;
                    Float Nnew = p.N + gamma * p.M;
                    Float Rnew = p.radius * std::sqrt(Nnew / (p.N + p.M));
                    Spectrum Phi;
                    for (int j = 0; j < Spectrum::nSamples; ++j)
                        Phi[j] = p.Phi[j];
                    p.tau = (p.tau + p.vp.beta * Phi) * (Rnew * Rnew) /
                            (p.radius * p.radius);
                    p.N = Nnew;
                    p.radius = Rnew;
                    p.M = 0;
                    for (int j = 0; j < Spectrum::nSamples; ++j)
                        p.Phi[j] = (Float) 0;
                }
                // Reset _VisiblePoint_ in pixel
                p.vp.beta = 0.0f;
                p.vp.bsdf = nullptr;
                p.vp.phase = nullptr;
            }, nPixels, 4096);
        }

        // Periodically store SPPM image in film and write image
        if (iter + 1 == nIterations || ((iter + 1) % writeFrequency) == 0) {
            int x0 = pixelBounds.pMin.x;
            int x1 = pixelBounds.pMax.x;
            uint64_t Np = (uint64_t) (iter + 1) * (uint64_t) photonsPerIteration;
            std::unique_ptr<Spectrum[]> image(new Spectrum[pixelBounds.Area()]);
            int offset = 0;
            for (int y = pixelBounds.pMin.y; y < pixelBounds.pMax.y; ++y) {
                for (int x = x0; x < x1; ++x) {
                    // Compute radiance _L_ for SPPM pixel _pixel_
                    const SPPMPixel &pixel =
                            pixels[(y - pixelBounds.pMin.y) * (x1 - x0) + (x - x0)];

                    // Contribute only if visible point should be rendered or not
                    Spectrum L = pixel.Ld / (iter + 1);
                    L += pixel.tau / (Np * Pi * pixel.radius * pixel.radius);
                    image[offset++] = L;
                }
            }
            camera->film->SetImage(image.get());
            camera->film->WriteImage();
            // Write SPPM radius image, if requested
            if (getenv("SPPM_RADIUS")) {
                std::unique_ptr<Float[]> rimg(
                        new Float[3 * pixelBounds.Area()]);
                Float minrad = 1e30f, maxrad = 0;
                for (int y = pixelBounds.pMin.y; y < pixelBounds.pMax.y; ++y) {
                    for (int x = x0; x < x1; ++x) {
                        const SPPMPixel &p =
                                pixels[(y - pixelBounds.pMin.y) * (x1 - x0) +
                                       (x - x0)];
                        minrad = std::min(minrad, p.radius);
                        maxrad = std::max(maxrad, p.radius);
                    }
                }
                fprintf(stderr,
                        "iterations: %d (%.2f s) radius range: %f - %f\n",
                        iter + 1, progress.ElapsedMS() / 1000., minrad, maxrad);
                int offset = 0;
                for (int y = pixelBounds.pMin.y; y < pixelBounds.pMax.y; ++y) {
                    for (int x = x0; x < x1; ++x) {
                        const SPPMPixel &p =
                                pixels[(y - pixelBounds.pMin.y) * (x1 - x0) +
                                       (x - x0)];
                        Float v = 1.f - (p.radius - minrad) / (maxrad - minrad);
                        rimg[offset++] = v;
                        rimg[offset++] = v;
                        rimg[offset++] = v;
                    }
                }
                Point2i res(pixelBounds.pMax.x - pixelBounds.pMin.x,
                            pixelBounds.pMax.y - pixelBounds.pMin.y);
                WriteImage("sppm_radius.png", rimg.get(), pixelBounds, res);
            }
        }
    }
    progress.Done();
}

Integrator *CreateVolSPPMIntegrator(const ParamSet &params,
                                    std::shared_ptr<const Camera> camera) {
    int nIterations =
            params.FindOneInt("iterations",
                              params.FindOneInt("numiterations", 64));
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int photonsPerIter = params.FindOneInt("photonsperiteration", -1);
    int writeFreq = params.FindOneInt("imagewritefrequency", 1 << 31);
    Float radius = params.FindOneFloat("radius", 1.f);
    if (PbrtOptions.quickRender) nIterations = std::max(1, nIterations / 16);

    bool renderSurfaces = params.FindOneBool("rendersurfaces", true);
    bool renderMedia = params.FindOneBool("rendermedia", true);

    return new VolSPPMIntegrator(camera, nIterations, photonsPerIter, maxDepth,
                                 radius, writeFreq,
                                 renderSurfaces, renderMedia);
}

}  // namespace pbrt
