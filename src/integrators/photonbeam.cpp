
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
#include "integrators/photonbeam.h"
#include "scene.h"
#include "imageio.h"
#include "paramset.h"
#include "progressreporter.h"
#include "sampling.h"
#include "samplers/halton.h"
#include "core/photonbeambvh.h"
#include <base/mutex.h> // from glog

namespace pbrt {

STAT_RATIO(
        "Progressive Photon Beams/Visible points checked per photon "
                "intersection",
        visiblePointsChecked, totalPhotonSurfaceInteractions);
STAT_COUNTER("Progressive Photon Beams/Photon paths followed",
             photonPaths);
STAT_COUNTER("Progressive Photon Beams/Total photon medium interactions",
             totalPhotonMediumInteractions);
STAT_COUNTER("Progressive Photon Beams/Total visible points in participating media",
             totalVisibleMediumPoints);
STAT_COUNTER("Progressive Photon Beams/Total visible points on surfaces",
             totalVisibleSurfacePoints);
STAT_INT_DISTRIBUTION(
        "Progressive Photon Beams/Grid cells per visible point",
        gridCellsPerVisiblePoint);
STAT_MEMORY_COUNTER("Memory/SPPM Pixels", pixelMemoryBytes);
STAT_FLOAT_DISTRIBUTION("Memory/SPPM BSDF and Grid Memory", memoryArenaMB);

// SPPM Local Definitions
struct PhotonBeamPixel {
    // PhotonBeamPixel Public Methods
    PhotonBeamPixel() {}

    // PhotonBeamPixel Public Data
    Spectrum Ld;
};

struct PhotonBeamPixelListNode {
    PhotonBeamPixel *pixel;
    PhotonBeamPixelListNode *next;
};

Float Determinant(const Vector3f &a, const Vector3f &b, const Vector3f &c) {
    // a b c
    // d e f
    // g h i
    return a.x * b.y * c.z + a.y * b.z * c.x + a.z * b.x * c.y
        - (a.z * b.y * c.x + a.y * b.x * c.z + a.x * b.z * c.y);
}

bool ComputeClosestPoints(const Point3f &a0, const Point3f &a1,
                          const Point3f &b0, const Point3f &b1,
                          Point3f &aClosest, Point3f &bClosest) {
    Vector3f A = a1 - a0;
    Vector3f B = b1 - b0;
    Float magA = A.Length();
    Float magB = B.Length();

    if (magA == 0.0f) {
        aClosest = a0;
        if (magB == 0.0f) {
            bClosest = b0;
            return true;
        }
        // Line segment A is a point, project A on B and clamp to end points
        B /= magB;

        A = a0 - b0;
        Float dot = Dot(A, B);
        bClosest = b0 + B * Clamp(dot, 0.0f, magB);

        return true;
    } else if (magB == 0.0f) {
        bClosest = b0;
        // Line segment B is a point, project B on A and clamp to end points
        A /= magA;

        B = b0 - a0;
        Float dot = Dot(A, B);
        aClosest = a0 + A * Clamp(dot, 0.0f, magA);

        return true;
    }

    A /= magA;
    B /= magB;

    const Vector3f cross = Cross(A, B);
    const Float denom = cross.LengthSquared();

    // If lines are parallel (denom=0) test if lines overlap.
    // If they don't overlap then there is a closest point solution.
    // If they do overlap, there are infinite closest positions,
    // but there is a closest distance
    if (denom == Float(0.0f)) {
        Float d0 = Dot(A, b0 - a0);
        Float d1 = Dot(A, b1 - a0);

        //# Is segment B before A?
        if (d0 <= 0 && d1 <= 0) {
            if (fabs(d0) < fabs(d1)) {
                aClosest = a0;
                bClosest = b1;
            }
            aClosest = a0;
            bClosest = b1;
        }

        // Is segment B after A?
        else if (d0 >= magA && d1 >= magA) {
            if (fabs(d0) < fabs(d1)) {
                aClosest = a1;
                bClosest = b1;
            }
            aClosest = a1;
            bClosest = b1;
        }

        return false;
    }

    // Lines criss-cross: Calculate the projected closest points
    const Vector3f t = b0 - a0;
    const Float detA = Determinant(t, B, cross);
    const Float detB = Determinant(t, A, cross);

    const Float t0 = detA / denom;
    const Float t1 = detB / denom;

    Point3f pA = a0 + A * t0;
    Point3f pB = b0 + B * t1;

    if (t0 < 0) pA = a0;
    else if (t0 > magA) pA = a1;

    // clamp projection of A
    if (t0 < 0 || t0 > magA) {
        Float dot = Clamp(Dot(B, pA - b0), 0.0f, magB);
        pB = b0 + B * dot;
    }

    if (t1 < 0 || t1 > magB) {
        Float dot = Clamp(Dot(A, pB - a0), 0.0f, magA);
        pA = a0 + A * dot;
    }

    aClosest = pA;
    bClosest = pB;
    return true;
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

    void Reset(uint64_t haltonIndex) {
        HaltonIndex = haltonIndex;
        haltonDim = 0;
        rng = RNG(haltonIndex);
    }

    std::unique_ptr<Sampler> Clone(int seed) override {
        return nullptr;
    }

private:
    uint64_t HaltonIndex;
    uint32_t haltonDim;
    RNG rng;
};

void TracePhotonBeamRecursive(Ray photonRay, int depth, Spectrum beta,
                              Sampler &localSampler, const Scene &scene,
                              const int MaxDepth, MemoryArena &arena, const Float BeamRadius,
                              std::vector<std::shared_ptr<PhotonBeam>> &localPhotonBeams) {
    SurfaceInteraction isect;
    for (; depth < MaxDepth; ++depth) {
        if (!scene.Intersect(photonRay, &isect)) break;

        Spectrum betaStart;

        MediumInteraction mi;
        Spectrum betaMedium(1.0f);
        if (photonRay.medium) betaMedium = photonRay.medium->Sample(photonRay, localSampler, arena, &mi);
        if (beta.IsBlack()) break;

        // Handle an interaction with a medium or a surface
        if (mi.IsValid()) {
            ++totalPhotonMediumInteractions;

            // Sample new photon ray direction from volume
            Vector3f wo = -photonRay.d, wi;
            mi.phase->Sample_p(wo, &wi, localSampler.Get2D());
            Ray mediumScatteredRay = mi.SpawnRay(wi);

            Spectrum scatteredBeta = beta * photonRay.medium->Tr(photonRay, localSampler);
            TracePhotonBeamRecursive(mediumScatteredRay, depth + 1, scatteredBeta,
                                     localSampler, scene,
                                     MaxDepth, arena, BeamRadius,
                                     localPhotonBeams);

        }

        ++totalPhotonSurfaceInteractions;
        if (photonRay.medium) betaMedium = photonRay.medium->Tr(photonRay, localSampler);
        auto beam = std::make_shared<PhotonBeam>(photonRay.o, isect.p, BeamRadius,
                                                 betaStart, betaMedium * beta);
        localPhotonBeams.push_back(beam);
        // Sample new photon ray direction from surface

        // Compute BSDF at photon intersection point
        isect.ComputeScatteringFunctions(photonRay, arena, true,
                                         TransportMode::Importance);
        if (!isect.bsdf) {
            --depth;
            photonRay = isect.SpawnRay(photonRay.d);
            continue;
        }
        const BSDF &isectBSDF = *isect.bsdf;

        // Sample BSDF _fr_ and direction _wi_ for reflected photon
        Vector3f wo = -photonRay.d, wi;
        Float pdf;
        BxDFType flags;

        // Generate _bsdfSample_ for outgoing photon sample
        Spectrum fr = isectBSDF.Sample_f(wo, &wi, localSampler.Get2D(), &pdf,
                                          BSDF_ALL, &flags);
        if (fr.IsBlack() || pdf == 0.f) break;
        Spectrum betaNew = betaMedium * beta * fr * AbsDot(wi, isect.shading.n) / pdf;

        photonRay = (RayDifferential) isect.SpawnRay(wi);

        // Possibly terminate photon path with Russian roulette
        Float q = std::max((Float) 0, 1 - betaNew.y() / beta.y());
        if (localSampler.Get1D() < q) break;
        beta = betaNew / (1 - q);
    }
}

// SPPM Method Definitions
void PhotonBeamIntegrator::Render(const Scene &scene) {
    ProfilePhase p(Prof::IntegratorRender);
    // Initialize _pixelBounds_ and _pixels_ array for SPPM
    Bounds2i pixelBounds = camera->film->croppedPixelBounds;
    int nPixels = pixelBounds.Area();
    std::unique_ptr<PhotonBeamPixel[]> pixels(new PhotonBeamPixel[nPixels]);
    const Float invSqrtSPP = 1.f / std::sqrt(nIterations);
    pixelMemoryBytes = nPixels * sizeof(PhotonBeamPixel);
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

    Float currentBeamRadius = initialBeamRadius;
    for (int iter = 0; iter < nIterations; ++iter) {
        // Counter to extend HaltonSampler values beyond 1000
        uint32_t iterNumPixels = 0;

        std::vector<std::shared_ptr<PhotonBeam>> photonBeams;

        // Trace photons and accumulate contributions
        {
            ProfilePhase _(Prof::SPPMPhotonPass);
            std::vector<MemoryArena> photonShootArenas(MaxThreadIndex());
            const int NumTasks = 1; //NumSystemCores();

            const int ChunkSize = photonsPerIteration / NumTasks;
            const int Remainder = photonsPerIteration - ChunkSize * NumTasks;
            std::vector<std::vector<std::shared_ptr<PhotonBeam>>> localPhotonBeamVector(NumTasks);

            ParallelFor([&](const uint32_t TaskIndex) {
                MemoryArena &arena = photonShootArenas[ThreadIndex];

                std::vector<std::shared_ptr<PhotonBeam>> &localPhotonBeams = localPhotonBeamVector[TaskIndex];
                const uint32_t PhotonStartIndex = ChunkSize * TaskIndex;

                int localChunkSize = ChunkSize;
                if (TaskIndex + 1 == NumTasks) localChunkSize += Remainder;

                for (uint32_t photonIndex = PhotonStartIndex;
                     photonIndex < PhotonStartIndex + localChunkSize; ++photonIndex) {
                    // Follow photon path for _photonIndex_
                    const uint64_t HaltonIndex =
                            (uint64_t) iter * (uint64_t) photonsPerIteration +
                            photonIndex;
                    AwesomeHaltonSampler localSampler(HaltonIndex+1);
                    std::vector<Float> samples;
                    for (uint32_t i = 0; i < 1000; ++i) samples.push_back(localSampler.Get1D());

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
                    if (pdfPos == 0 || pdfDir == 0 || Le.IsBlack()) continue;
                    Spectrum beta = (AbsDot(nLight, photonRay.d) * Le) /
                                    (lightPdf * pdfPos * pdfDir);
                    if (beta.IsBlack()) continue;
                    TracePhotonBeamRecursive(photonRay, 0, beta,
                                            localSampler, scene,
                                            maxDepth, arena, currentBeamRadius,
                                            localPhotonBeams);
                    arena.Reset();
                }
            }, NumTasks, 1);
            progress.Update();
            photonPaths += photonsPerIteration;

            // Merge local photon beams into global
            long size = 0;
            for (auto &localPhotonBeams : localPhotonBeamVector)
                size += localPhotonBeams.size();
            photonBeams.reserve(size);
            for (auto &localPhotonBeams : localPhotonBeamVector) {
                for (uint32_t i = 0; i < localPhotonBeams.size(); ++i)
                    photonBeams.push_back(localPhotonBeams[i]);
                localPhotonBeams.erase(localPhotonBeams.begin(),
                                       localPhotonBeams.end());
            }
        }
        PhotonBeamBVH photonBeamBVH(std::move(photonBeams));

        // Cast camera rays and gather contributions from photon beams
        std::vector<MemoryArena> perThreadArenas(MaxThreadIndex());
        {
            ProfilePhase _(Prof::SPPMCameraPass);
            ParallelFor2D([&](Point2i tile) {
                MemoryArena &arena = perThreadArenas[ThreadIndex];
                // Follow camera paths for _tile_ in im age for VolSPPM
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

                    // Get _PhotonBeamPixel_ for _pPixel_
                    Point2i pPixelO = Point2i(pPixel - pixelBounds.pMin);
                    int pixelOffset =
                            pPixelO.x +
                            pPixelO.y * (pixelBounds.pMax.x - pixelBounds.pMin.x);
                    PhotonBeamPixel &pixel = pixels[pixelOffset];
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

                        ++totalPhotonSurfaceInteractions;

                        Spectrum mediumBeta(1.0f);
                        if (ray.medium) mediumBeta = ray.medium->Tr(ray, localSampler);

                        if (renderMedia) {
                            std::vector<std::shared_ptr<PhotonBeam>> beams = photonBeamBVH.Intersect(ray);
                            for (std::shared_ptr<PhotonBeam> const& beam : beams) {
                                // Add contribution of photon beam along camera ray
                                Point3f rayClose, beamClose;
                                if (ComputeClosestPoints(ray.o, isect.p, beam->start, beam->end, rayClose, beamClose)) {
                                    const Float MaxDistance = currentBeamRadius + beam->radius;
                                    Float distance = (rayClose - beamClose).Length();
                                    if (distance < MaxDistance) {
                                        Float r = distance / MaxDistance;
                                        pixel.Ld += 1e-5 * beam->powerEnd * sqrt(1.0f - r * r);
                                    }
                                }
                            }
                        }

                        beta *= mediumBeta;

                        // Process SPPM camera ray surface intersection
                        // Compute BSDF at SPPM camera ray intersection
                        isect.ComputeScatteringFunctions(ray, arena, true);
                        if (!isect.bsdf) {
                            ray = isect.SpawnRay(ray.d);
                            --depth;
                        } else {
                            if (!renderSurfaces) {
                                break;
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
                                ray = (RayDifferential) isect.SpawnRay(wi);
                            }
                        }

                        if (beta.y() < 0.25) {
                            Float continueProb =
                                    std::min((Float) 1, beta.y());
                            if (localSampler.Get1D() > continueProb) break;
                            beta /= continueProb;
                        }
                    }
                }
            }, nTiles);
        }
        progress.Update();

        // update beam radius according to variance/bias rule
        currentBeamRadius = currentBeamRadius * (Float(iter + alpha) / Float(iter + 1));

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
                    const PhotonBeamPixel &pixel =
                            pixels[(y - pixelBounds.pMin.y) * (x1 - x0) + (x - x0)];

                    // Contribute only if visible point should be rendered or not
                    Spectrum L = pixel.Ld / (iter + 1);
                    image[offset++] = L;
                }
            }
            camera->film->SetImage(image.get());
            camera->film->WriteImage();
        }
    }
    progress.Done();
}

Integrator *CreatePhotonBeamIntegrator(const ParamSet &params,
                                    std::shared_ptr<const Camera> camera) {
    int nIterations =
            params.FindOneInt("iterations",
                              params.FindOneInt("numiterations", 64));
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int photonsPerIter = params.FindOneInt("photonsperiteration", -1);
    int writeFreq = params.FindOneInt("imagewritefrequency", 1 << 31);
    Float radius = params.FindOneFloat("initialbeamradius", 1.f);
    Float alpha = params.FindOneFloat("alpha", 0.5f);
    if (PbrtOptions.quickRender) nIterations = std::max(1, nIterations / 16);

    bool renderSurfaces = params.FindOneBool("rendersurfaces", true);
    bool renderMedia = params.FindOneBool("rendermedia", true);

    return new PhotonBeamIntegrator(camera, nIterations, photonsPerIter, maxDepth,
                                    radius, alpha,
                                    writeFreq,
                                    renderSurfaces, renderMedia);
}

}  // namespace pbrt
