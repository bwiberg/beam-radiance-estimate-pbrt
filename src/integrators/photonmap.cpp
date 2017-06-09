
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

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


// integrators/photonmap.cpp*
#include <memory>
#include <array>
#include <samplers/random.h>
#include "base/mutex.h"
#include "samplers/halton.h"
#include "integrators/photonmap.h"
#include "scene.h"
#include "progressreporter.h"
#include "paramset.h"
#include "camera.h"

namespace pbrt {
STAT_COUNTER("Photon Mapping/Total Photon Surface Interactions ", totalPhotonSurfaceInteractions);
STAT_COUNTER("Photon Mapping/Total Photon Surface Interactions ", totalPhotonVolumeInteractions);
STAT_COUNTER("Photon Mapping/Photon paths followed", photonPaths);
STAT_COUNTER("Photon Mapping/Total indirect photons", indirectPhotonCounter);
STAT_COUNTER("Photon Mapping/Total caustic photons", causticPhotonCounter);
STAT_COUNTER("Photon Mapping/Total radiance photons", radiancePhotonCounter);
STAT_COUNTER("Photon Mapping/Total volume photons", volumePhotonCounter);
STAT_MEMORY_COUNTER("Memory/Photon maps", photonMapBytes);

template<size_t N>
std::array<Float, N> GetSamples1D(Sampler &sampler);

template<size_t N>
std::array<Point2f, N> GetSamples2D(Sampler &sampler);

struct Photon {
    Photon(const Point3f &pp, const Spectrum &wt, const Vector3f &w)
            : p(pp), alpha(wt), wi(w) {}

    Photon() {}

    Point3f p;
    Spectrum alpha;
    Vector3f wi;
};

struct RadiancePhoton {
    RadiancePhoton(const Point3f &pp, const Normal3f &nn)
            : p(pp), n(nn), Lo(0.f) {}

    RadiancePhoton() {}

    Point3f p;
    Normal3f n;
    Spectrum Lo;
};


struct PhotonProcess {
    // PhotonProcess Public Methods
    PhotonProcess(uint32_t mp, ClosePhoton *buf);

    void operator()(const Point3f &p, const Photon &photon, Float dist2,
                    Float &maxDistSquared);

    ClosePhoton *photons;
    uint32_t nLookup, nFound;
};


struct ClosePhoton {
    // ClosePhoton Public Methods
    ClosePhoton(const Photon *p = NULL, Float md2 = INFINITY)
            : photon(p), distanceSquared(md2) {}

    bool operator<(const ClosePhoton &p2) const {
        return distanceSquared == p2.distanceSquared ?
               (photon < p2.photon) : (distanceSquared < p2.distanceSquared);
    }

    const Photon *photon;
    Float distanceSquared;
};


PhotonProcess::PhotonProcess(uint32_t mp, ClosePhoton *buf)
        : photons(buf), nLookup(mp), nFound(0) {}


struct RadiancePhotonProcess {
    // RadiancePhotonProcess Methods
    RadiancePhotonProcess(const Normal3f &nn)
            : n(nn) {
        photon = NULL;
    }

    void operator()(const Point3f &p, const RadiancePhoton &rp,
                    Float distSquared, Float &maxDistSquared) {
        if (Dot(rp.n, n) > 0) {
            photon = &rp;
            maxDistSquared = distSquared;
        }
    }

    const Normal3f &n;
    const RadiancePhoton *photon;
};


inline Float kernel(const Photon *photon, const Point3f &p, Float maxDist2);

static Spectrum LPhoton(KdTree<Photon> const *map, int nPaths, int nLookup,
                        ClosePhoton *lookupBuf, BSDF *bsdf, RNG &rng, const SurfaceInteraction &isect,
                        const Vector3f &w, Float maxDistSquared);

static Spectrum EPhoton(KdTree<Photon> const *map, int count, int nLookup,
                        ClosePhoton *lookupBuf, Float maxDist2, const Point3f &p, const Normal3f &n);

// PhotonIntegrator Local Deff
// ifinitions
inline bool unsuccessful(uint32_t needed, uint32_t found, uint32_t shot) {
    return (found < needed && (found == 0 || found < shot / 1024));
}


inline void PhotonProcess::operator()(const Point3f &p,
                                      const Photon &photon, Float distSquared, Float &maxDistSquared) {
    if (nFound < nLookup) {
        // Add photon to unordered array of photons
        photons[nFound++] = ClosePhoton(&photon, distSquared);
        if (nFound == nLookup) {
            std::make_heap(&photons[0], &photons[nLookup]);
            maxDistSquared = photons[0].distanceSquared;
        }
    } else {
        // Remove most distant photon from heap and add new photon
        std::pop_heap(&photons[0], &photons[nLookup]);
        photons[nLookup - 1] = ClosePhoton(&photon, distSquared);
        std::push_heap(&photons[0], &photons[nLookup]);
        maxDistSquared = photons[0].distanceSquared;
    }
}


inline Float kernel(const Photon *photon, const Point3f &p,
                    Float maxDist2) {
    Float s = (1.f - DistanceSquared(photon->p, p) / maxDist2);
    return 3.f * InvPi * s * s;
}


Spectrum LPhoton(KdTree<Photon> const *map, int nPaths, int nLookup,
                 ClosePhoton *lookupBuf, BSDF *bsdf, Sampler &sampler,
                 const SurfaceInteraction &isect, const Vector3f &wo, Float maxDist2) {
    Spectrum L(0.);
    BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
                                    BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
    if (map && bsdf->NumComponents(nonSpecular) > 0) {
        // Do photon map lookup at intersection Point3f
        PhotonProcess proc(nLookup, lookupBuf);
        map->Lookup(isect.p, proc, maxDist2);

        // Estimate reflected radiance due to incident photons
        ClosePhoton *photons = proc.photons;
        int nFound = proc.nFound;
        Normal3f Nf = Faceforward(isect.n, wo);
        if (bsdf->NumComponents(BxDFType(BSDF_REFLECTION |
                                         BSDF_TRANSMISSION | BSDF_GLOSSY)) > 0) {
            // Compute exitant radiance from photons for glossy surface
            for (int i = 0; i < nFound; ++i) {
                const Photon *p = photons[i].photon;
                Float k = kernel(p, isect.p, maxDist2);
                L += (k / (nPaths * maxDist2)) * bsdf->f(wo, p->wi) *
                     p->alpha;
            }
        } else {
            // Compute exitant radiance from photons for diffuse surface
            Spectrum Lr(0.), Lt(0.);
            for (int i = 0; i < nFound; ++i) {
                if (Dot(Nf, photons[i].photon->wi) > 0.f) {
                    Float k = kernel(photons[i].photon, isect.p, maxDist2);
                    Lr += (k / (nPaths * maxDist2)) * photons[i].photon->alpha;
                } else {
                    Float k = kernel(photons[i].photon, isect.p, maxDist2);
                    Lt += (k / (nPaths * maxDist2)) * photons[i].photon->alpha;
                }
            }

            // TODO VALIDATE
            constexpr int sqrtSamples = 6;
            constexpr int nSamples = sqrtSamples * sqrtSamples;
            auto s1 = GetSamples2D<nSamples>(sampler);
            auto s2 = GetSamples2D<nSamples>(sampler);

            L += Lr * bsdf->rho(wo, nSamples, s1.data(), BxDFType(BSDF_ALL_TYPES | BSDF_REFLECTION)) * InvPi +
                 Lt * bsdf->rho(wo, nSamples, s2.data(), BxDFType(BSDF_ALL_TYPES | BSDF_TRANSMISSION)) * InvPi;
        }
    }
    return L;
}


Spectrum EPhoton(KdTree<Photon> const *map, int count, int nLookup,
                 ClosePhoton *lookupBuf, Float maxDist2, const Point3f &p,
                 const Normal3f &n) {
    if (!map) return 0.f;
    // Lookup nearby photons at irradiance computation Point3f
    PhotonProcess proc(nLookup, lookupBuf);
    Float md2 = maxDist2;
    map->Lookup(p, proc, md2);
    DCHECK(md2 > 0.f);

    // Accumulate irradiance value from nearby photons
    if (proc.nFound == 0) return Spectrum(0.f);
    ClosePhoton *photons = proc.photons;
    Spectrum E(0.);
    for (uint32_t i = 0; i < proc.nFound; ++i)
        if (Dot(n, photons[i].photon->wi) > 0.)
            E += photons[i].photon->alpha;
    Float denominator = (count * md2 * M_PI);
    if (denominator == 0.f)
        return Spectrum(0.f);
    return E / denominator;
}


// PhotonIntegrator Method Definitions
PhotonIntegrator::PhotonIntegrator(std::shared_ptr<const Camera> camera,
                                   std::shared_ptr<Sampler> sampler,
                                   const Bounds2i &pixelBounds,
                                   int ncaus, int nind, int nvolume, bool reqphotons,
                                   int nl, int mdepth, int mphodepth, float mdist, bool fg,
                                   int gs, float ga) :
        SamplerIntegrator(camera, sampler, pixelBounds),
        nCausticPhotonsWanted(ncaus),
        nIndirectPhotonsWanted(nind),
        nVolumePhotonsWanted(nvolume),
        nLookup(nl),
        requirePhotons(reqphotons),
        maxSpecularDepth(mdepth),
        maxPhotonDepth(mphodepth),
        maxDistSquared(mdist * mdist),
        finalGather(fg),
        cosGatherAngle(cos(Radians(ga))),
        gatherSamples(gs),
        nCausticPaths(nIndirectPaths = 0),
        causticMap(nullptr),
        indirectMap(nullptr),
        radianceMap(nullptr) {}


PhotonIntegrator::~PhotonIntegrator() {
    delete causticMap;
    delete indirectMap;
    delete radianceMap;
}


void PhotonIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    int const numLights = scene.lights.size();
    if (numLights == 0) return;

    ComputeLightSamples(scene, sampler);
    ComputeLightPhotonDistrib(scene);

    // Declare shared variables for photon shooting
    int nDirectPaths = 0;
    uint32_t nshot = 0;
    std::vector<Photon> causticPhotons, directPhotons, indirectPhotons, volumePhotons;
    std::vector<RadiancePhoton> radiancePhotons;
    causticPhotons.reserve(nCausticPhotonsWanted);
    indirectPhotons.reserve(nIndirectPhotonsWanted);
    volumePhotons.reserve(nVolumePhotonsWanted);
    std::vector<Spectrum> rpReflectances, rpTransmittances;

    ShootPhotons(scene,
                 causticPhotons, directPhotons, indirectPhotons, volumePhotons,
                 radiancePhotons, rpReflectances, rpTransmittances,
                 nDirectPaths, nshot);

    BuildPhotonMaps(directPhotons, causticPhotons, indirectPhotons, volumePhotons);
//
//    if (finalGather && radiancePhotons.size() > 0) {
//        ComputePhotonRadiances(radiancePhotons,
//                               rpReflectances, rpTransmittances,
//                               directMap,
//                               nDirectPaths);
//    }
//
//    {
//        // Update PhotonMapping stats
//        indirectPhotonCounter = static_cast<int64_t>(indirectPhotons.size());
//        causticPhotonCounter = static_cast<int64_t>(causticPhotons.size());
//        radiancePhotonCounter = static_cast<int64_t>(radiancePhotons.size());
//        volumePhotonCounter = 0;
//        photonMapBytes = sizeof(Photon) * (indirectPhotonCounter + causticPhotonCounter)
//                         + sizeof(RadiancePhoton) * radiancePhotonCounter;
//    }
//
//    delete directMap;
}


Spectrum PhotonIntegrator::Li(const RayDifferential &ray,
                              const Scene &scene,
                              Sampler &sampler,
                              MemoryArena &arena,
                              int depth) const {
    Spectrum L(0.0f);

    // Find closest ray intersection or return background radiance
    SurfaceInteraction isect;
    if (!scene.Intersect(ray, &isect)) {
        for (const auto &light : scene.lights) L += light->Le(ray);
        return L;
    }

    const Point3f p = isect.p;

    const uint32_t nIndirSamplePhotons = 50;
    PhotonProcess proc(nIndirSamplePhotons,
                       arena.Alloc<ClosePhoton>(nIndirSamplePhotons));
    Float md2 = 10 * maxDistSquared;
    if (indirectMap)
        indirectMap->Lookup(p, proc, md2);
    if (proc.nFound > 0)
        L = Spectrum(proc.nFound * 1.0f);

    return L;
//    ProfilePhase prof(Prof::SamplerIntegratorLi);
//
//    Spectrum L(0.0f);
//
//    // Find closest ray intersection or return background radiance
//    SurfaceInteraction isect;
//    if (!scene.Intersect(ray, &isect)) {
//        for (const auto &light : scene.lights) L += light->Le(ray);
//        return L;
//    }
//
//    // Compute scattering functions for surface interaction
//    isect.ComputeScatteringFunctions(ray, arena);
//    if (!isect.bsdf)
//        return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
//
//    Vector3f wo = -ray.d;
//
//    // Compute emitted light if ray hit an area light source
//    L += isect.Le(wo);
//
//    // Evaluate BSDF at hit Point3f
//    const Point3f &p = isect.p;
//    const Normal3f &n = isect.n;
//
//    bool handleMedia = false;
//    L += UniformSampleAllLights(isect, scene, arena, sampler, nLightSamples, handleMedia);
//    // Compute caustic lighting for photon map integrator
//    ClosePhoton *lookupBuf = arena.Alloc<ClosePhoton>(nLookup);
//    L += LPhoton(causticMap, nCausticPaths, nLookup, lookupBuf, isect.bsdf,
//                 sampler, isect, wo, maxDistSquared);
//
//    // Compute indirect lighting for photon map integrator
//    if (finalGather && indirectMap != NULL) {
//#if 1
//        // Do one-bounce final gather for photon map
//        BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
//                                        BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
//        if (isect.bsdf->NumComponents(nonSpecular) > 0) {
//            // Find indirect photons around Point3f  for importance sampling
//            const uint32_t nIndirSamplePhotons = 50;
//            PhotonProcess proc(nIndirSamplePhotons,
//                               arena.Alloc<ClosePhoton>(nIndirSamplePhotons));
//            Float searchDist2 = maxDistSquared;
//            while (proc.nFound < nIndirSamplePhotons) {
//                Float md2 = searchDist2;
//                proc.nFound = 0;
//                indirectMap->Lookup(p, proc, md2);
//                searchDist2 *= 2.f;
//            }
//
//            // Copy photon directions to local array
//            Vector3f *photonDirs = arena.Alloc<Vector3f>(nIndirSamplePhotons);
//            for (uint32_t i = 0; i < nIndirSamplePhotons; ++i)
//                photonDirs[i] = proc.photons[i].photon->wi;
//
//            // Use BSDF to do final gathering
//            Spectrum Li = 0.;
//            for (int i = 0; i < gatherSamples; ++i) {
//                // Sample random direction from BSDF for final gather ray
//                Vector3f wi;
//                Float pdf;
//                Point2f u_bsdf;
//                Spectrum fr = isect.bsdf->Sample_f(wo, &wi, u_bsdf,
//                                                   &pdf, BxDFType(BSDF_ALL & ~BSDF_SPECULAR));
//                if (fr.IsBlack() || pdf == 0.f) continue;
//                DCHECK(pdf >= 0.f);
//
//                // Trace BSDF final gather ray and accumulate radiance
//                Ray bounceRay = isect.SpawnRay(wi);
//                SurfaceInteraction gatherIsect;
//                if (scene.Intersect(bounceRay, &gatherIsect)) {
//                    // Compute exitant radiance _Lindir_ using radiance photons
//                    Spectrum Lindir = 0.f;
//                    Normal3f nGather = gatherIsect.n;
//                    nGather = Faceforward(nGather, -bounceRay.d);
//                    RadiancePhotonProcess proc(nGather);
//                    Float md2 = INFINITY;
//                    radianceMap->Lookup(gatherIsect.p, proc, md2);
//                    if (proc.photon != NULL)
//                        Lindir = proc.photon->Lo;
//
//                    // TODO verify
//                    if (bounceRay.medium)
//                        Lindir *= bounceRay.medium->Tr(bounceRay, sampler);
//
//
//                    // Compute MIS weight for BSDF-sampled gather ray
//
//                    // Compute PDF for photon-sampling of direction _wi_
//                    Float photonPdf = 0.f;
//                    Float conePdf = UniformConePdf(cosGatherAngle);
//                    for (uint32_t j = 0; j < nIndirSamplePhotons; ++j)
//                        if (Dot(photonDirs[j], wi) > .999f * cosGatherAngle)
//                            photonPdf += conePdf;
//                    photonPdf /= nIndirSamplePhotons;
//                    Float wt = PowerHeuristic(gatherSamples, pdf, gatherSamples, photonPdf);
//                    Li += fr * Lindir * (AbsDot(wi, n) * wt / pdf);
//                }
//            }
//            L += Li / gatherSamples;
//
//            // Use nearby photons to do final gathering
//            Li = 0.;
//            for (int i = 0; i < gatherSamples; ++i) {
//                // Sample random direction using photons for final gather ray
//                // TODO
//                // TODO
//                // TODO
//                // TODO
//                // TODO what is BSDFSample equivalent?
//                // TODO => random sampler.2D
//                // TODO
//                // TODO
//                // TODO
//                Point2f u_cone = sampler.Get2D();
//                int photonNum = std::min((int) nIndirSamplePhotons - 1,
//                                         static_cast<int>(std::floor(sampler.Get1D() * nIndirSamplePhotons)));
//
//                // Sample gather ray direction from _photonNum_
//                Vector3f vx, vy;
//                CoordinateSystem(photonDirs[photonNum], &vx, &vy);
//                Vector3f wi = UniformSampleCone(u_cone, cosGatherAngle, vx, vy, photonDirs[photonNum]);
//                // Trace photon-sampled final gather ray and accumulate radiance
//                Spectrum fr = isect.bsdf->f(wo, wi);
//                if (fr.IsBlack()) continue;
//                RayDifferential bounceRay(p, wi);
//                SurfaceInteraction gatherIsect;
//                if (scene.Intersect(bounceRay, &gatherIsect)) {
//                    // Compute exitant radiance _Lindir_ using radiance photons
//                    Spectrum Lindir = 0.f;
//                    Normal3f nGather = gatherIsect.n;
//                    nGather = Faceforward(nGather, -bounceRay.d);
//                    RadiancePhotonProcess proc(nGather);
//                    Float md2 = INFINITY;
//                    radianceMap->Lookup(gatherIsect.p, proc, md2);
//                    if (proc.photon != NULL)
//                        Lindir = proc.photon->Lo;
//
//                    // TODO verify
//                    if (bounceRay.medium)
//                        Lindir *= bounceRay.medium->Tr(bounceRay, sampler);
//
//                    // Compute PDF for photon-sampling of direction _wi_
//                    Float photonPdf = 0.f;
//                    Float conePdf = UniformConePdf(cosGatherAngle);
//                    for (uint32_t j = 0; j < nIndirSamplePhotons; ++j)
//                        if (Dot(photonDirs[j], wi) > .999f * cosGatherAngle)
//                            photonPdf += conePdf;
//                    photonPdf /= nIndirSamplePhotons;
//
//                    // Compute MIS weight for photon-sampled gather ray
//                    Float bsdfPdf = isect.bsdf->Pdf(wo, wi);
//                    Float wt = PowerHeuristic(gatherSamples, photonPdf, gatherSamples, bsdfPdf);
//                    Li += fr * Lindir * AbsDot(wi, n) * wt / photonPdf;
//                }
//            }
//            L += Li / gatherSamples;
//        }
//#else
//        // for debugging / examples: use the photon map directly
//        Normal3f nn = Faceforward(n, -ray.d);
//        RadiancePhotonProcess proc(nn);
//        Float md2 = INFINITY;
//        radianceMap->Lookup(p, proc, md2);
//        if (proc.photon)
//            L += proc.photon->Lo;
//#endif
//    } else
//        L += LPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
//                     isect.bsdf, sampler, isect, wo, maxDistSquared);
//    if (depth + 1 < maxSpecularDepth) {
//        Vector3f wi;
//
//        // Trace rays for specular reflection and refraction
//        L += SpecularReflect(ray, isect, scene, sampler, arena, depth);
//        L += SpecularTransmit(ray, isect, scene, sampler, arena, depth);
//    }
//    return L;
}


void PhotonIntegrator::ComputeLightSamples(const Scene &scene, Sampler &sampler) {
    for (int i = 0; i < scene.lights.size(); ++i) {
        // Compute number of samples to use for each light
        nLightSamples.push_back(sampler.RoundCount(scene.lights[i]->nSamples));
    }
}


void PhotonIntegrator::ComputeLightPhotonDistrib(Scene const &scene) {
    // compute power-based probability distribution of lights
    auto const numLights = scene.lights.size();

    std::vector<Float> lightPbs(numLights);
    Float sumLightPowers = 0.0;
    for (int i = 0; i < numLights; ++i) {
        lightPbs[i] = scene.lights[i]->Power().MaxComponentValue();
        sumLightPowers += lightPbs[i];
    }
    DCHECK_GT(sumLightPowers, 0.0);
    for (int i = 0; i < numLights; ++i) {
        lightPbs[i] /= sumLightPowers;
    }
    lightDistr = std::unique_ptr<Distribution1D>(new Distribution1D(lightPbs.data(), numLights));
}


class AwesomeSampler : public Sampler {
public:
    AwesomeSampler(int64_t samplesPerPixel, uint64_t photonIndex, RNG &rng)
            : Sampler(samplesPerPixel),
              rng(rng),
              photonIndex(photonIndex),
              haltonDim(0) {}

    Float Get1D() override {
        if (haltonDim < 1000) {
            return RadicalInverse(haltonDim++, photonIndex);
        }
        return rng.UniformFloat();
    }

    Point2f Get2D() override {
        if (haltonDim + 1 <= 1000) {
            const int dim0 = haltonDim++;
            return Point2f(RadicalInverse(dim0, photonIndex),
                           RadicalInverse(haltonDim++, photonIndex));
        }
        return Point2f(rng.UniformFloat(), rng.UniformFloat());
    }

    std::unique_ptr<Sampler> Clone(int seed) override {
        return nullptr;
    }

private:
    RNG &rng;
    uint64_t photonIndex;
    int haltonDim;
};


void
PhotonIntegrator::ShootPhotons(Scene const &scene,
                               std::vector<Photon> &causticPhotons,
                               std::vector<Photon> &directPhotons,
                               std::vector<Photon> &indirectPhotons,
                               std::vector<Photon> &volumePhotons,
                               std::vector<RadiancePhoton> &radiancePhotons,
                               std::vector<Spectrum> &rpReflectances, std::vector<Spectrum> &rpTransmittances,
                               int &nDirectPaths, uint32_t &nShot) {
    ProgressReporter progress(nCausticPhotonsWanted + nIndirectPhotonsWanted + nVolumePhotonsWanted, "Shooting photons");

    Mutex mutex;
    int const nTasks = NumSystemCores();
    bool abortTasks = false;

    int64_t numProgress = 0;

    uint32_t totalPhotonsSent = 0;

    // Run parallel tasks for photon shooting
    ParallelFor([&](int i_task) {
        // Loop until desired number of photons met:
        // 1. shoot n photons
        // 2. append local direct, indirect and caustic photon vectors to global vectors

        // Declare local variables for _PhotonShootingTask_
        MemoryArena arena;
        std::vector<Photon> localDirectPhotons, localIndirectPhotons, localCausticPhotons, localVolumePhotons;
        std::vector<RadiancePhoton> localRadiancePhotons;
        uint32_t totalPaths = 0;
        bool causticDone = (nCausticPhotonsWanted == 0);
        bool indirectDone = (nIndirectPhotonsWanted == 0);
        bool volumeDone = (nVolumePhotonsWanted == 0);
        std::vector<Spectrum> localRpReflectances, localRpTransmittances;

        for (;;) {
            constexpr uint32_t BlockSize = 4096;

            uint32_t LocalPhotonStartIndex;
            {
                MutexLock lock(&mutex);
                LocalPhotonStartIndex = totalPhotonsSent;
                totalPhotonsSent += BlockSize;
            }
            // Follow photon paths for a block of samples
            for (uint32_t i = 0; i < BlockSize; ++i) {
                const uint32_t PhotonIndex = LocalPhotonStartIndex + i;
                RNG rng(PhotonIndex);
                AwesomeSampler sampler(0, PhotonIndex, rng);

                auto Get2Ds = [&](const size_t NumValues) {
                    Point2f *values = arena.Alloc<Point2f>(NumValues, false);
                    for (size_t j = 0; j < NumValues; ++j)
                        values[j] = sampler.Get2D();
                    return values;
                };

                Float const u0 = sampler.Get1D();
                Point2f const &u12 = sampler.Get2D();
                Point2f const &u34 = sampler.Get2D();
                //Stepsampler();

                // Choose light to shoot photon from
                Float lightPdf;
                int lightNum = lightDistr->SampleDiscrete(u0, &lightPdf);
                std::shared_ptr<Light> light = scene.lights[lightNum];

                // Generate _photonRay_ from light source and initialize _alpha_
                RayDifferential photonRay;
                Vector3f wi;
                Float pdfLightPos, pdfLightDir;
                Normal3f nLight;
                Spectrum Le = light->Sample_Le(u12, u34,
                                               (camera ? camera->shutterOpen : 0.0f), &photonRay, &nLight, &pdfLightPos,
                                               &pdfLightDir);

                if (pdfLightPos == 0.0f || pdfLightDir == 0.0f || Le.IsBlack()) continue;
                ++photonPaths;

                Spectrum alpha = (AbsDot(nLight, photonRay.d) * Le) / (pdfLightDir * pdfLightPos * lightPdf);
                if (!alpha.IsBlack()) {
                    // Follow photon path through scene and record intersections
                    bool specularPath = true;
                    SurfaceInteraction photonIsect;
                    int nIntersections = 0;
                    while (scene.Intersect(photonRay, &photonIsect)) {
                        ++nIntersections;

                        // Participating media interactions
                        MediumInteraction mi;
                        if (photonRay.medium) alpha *= photonRay.medium->Sample(photonRay, sampler, arena, &mi);
                        if (alpha.IsBlack()) break;

                        // Handle photon interaction with medium OR surface
                        photonIsect.ComputeScatteringFunctions(photonRay, arena);
                        BSDF *bsdf = photonIsect.bsdf;
                        if (!bsdf) {
                            photonRay = photonIsect.SpawnRay(photonRay.d);
                            nIntersections--;
                            continue;
                        }

                        Vector3f wo = -photonRay.d;
                        if (mi.IsValid()) {
                            ++totalPhotonVolumeInteractions;

                            Photon photon(mi.p, alpha, wo);
                            if (nIntersections > 1 && !volumeDone)
                                localVolumePhotons.push_back(photon);

                        } else {
                            ++totalPhotonSurfaceInteractions;

                            // Handle photon/surface interaction
                            constexpr BxDFType SpecularBxDFType = BxDFType(BSDF_REFLECTION |
                                                                           BSDF_TRANSMISSION |
                                                                           BSDF_SPECULAR);
                            const bool HasNonSpecular = bsdf->NumComponents() >
                                                        bsdf->NumComponents(SpecularBxDFType);
                            if (HasNonSpecular) {
                                // Deposit photon at surface
                                const Photon photon(photonIsect.p, alpha, wo);
                                bool hasDepositedPhoton = false;
                                if (specularPath && nIntersections > 1 && !causticDone) {
                                    hasDepositedPhoton = true;
                                    localCausticPhotons.push_back(photon);
                                } else {
                                    // Deposit either direct or indirect photon
                                    // stop depositing direct photons once indirectDone is true; don't
                                    // want to waste memory storing too many if we're going a long time
                                    // trying to get enough caustic photons desposited.
                                    if (nIntersections == 1 && !indirectDone && finalGather) {
                                        hasDepositedPhoton = true;
                                        localDirectPhotons.push_back(photon);
                                    } else if (nIntersections > 1 && !indirectDone) {
                                        hasDepositedPhoton = true;
                                        localIndirectPhotons.push_back(photon);
                                    }
                                }

                                // Possibly create radiance photon at photon intersection point
//                            if (hasDepositedPhoton && finalGather && sampler.Get1D() < .125f) {
//                                const Normal3f normal = Faceforward(photonIsect.n, -photonRay.d);
//                                localRadiancePhotons.push_back(RadiancePhoton(photonIsect.p, normal));
//
//                                constexpr int32_t SqrtSamples = 6;
//                                constexpr int32_t NumSamples = SqrtSamples * SqrtSamples;
//
//                                Point2f *s1 = Get2Ds(NumSamples);
//                                Point2f *s2 = Get2Ds(NumSamples);
//                                //std::array<Point2f, NumSamples> s1 = GetSamples2D<NumSamples>(*sampler);
//                                //std::array<Point2f, NumSamples> s2 = GetSamples2D<NumSamples>(*sampler);
//                                constexpr BxDFType ReflectanceBxDFType = BxDFType(BSDF_ALL_TYPES | BSDF_REFLECTION);
//                                localRadpRfls.push_back(bsdf->rho(NumSamples, s1, s2, ReflectanceBxDFType));
//
//                                s1 = Get2Ds(NumSamples);
//                                s2 = Get2Ds(NumSamples);
//                                constexpr BxDFType TransmittanceBxDFType = BxDFType(BSDF_ALL_TYPES | BSDF_TRANSMISSION);
//                                localRadpTrs.push_back(bsdf->rho(NumSamples, s1, s2, TransmittanceBxDFType));
//                            }
                            }
                        }
                        if (nIntersections >= this->maxPhotonDepth) break;

                        // Sample new photon ray direction
                        Float pdf;
                        BxDFType flags;
                        Point2f uphoton = sampler.Get2D();
                        Spectrum fr = bsdf->Sample_f(-photonRay.d, &wi, uphoton, &pdf, BSDF_ALL, &flags);
                        if (fr.IsBlack() || pdf == 0.f) break;
                        Spectrum anew = alpha * fr *
                                        AbsDot(wi, photonIsect.n) / pdfLightPos;

                        // Possibly terminate photon path with Russian roulette
                        Float continueProb = std::min(1.f, anew.y() / alpha.y());
                        Float sample = sampler.Get1D();
                        if (sample > continueProb)
                            break;
                        alpha = anew / continueProb;
                        specularPath &= ((flags & BSDF_SPECULAR) != 0);

                        if (indirectDone && !specularPath) break;
                        photonRay = photonIsect.SpawnRay(wi);
                    }
                }
                arena.Reset();
            }

            // Merge local photon data with data in _PhotonIntegrator_
            {
                MutexLock lock(&mutex);

                // If photons aren't required
                // => Give up if we're not storing enough photons
                // NOT guaranteed to ever complete if (requirePhotons == true)
                if (abortTasks)
                    return;
                if (!requirePhotons && nShot > 500000 &&
                    (unsuccessful(this->nCausticPhotonsWanted,
                                  causticPhotons.size(), BlockSize) ||
                     unsuccessful(this->nIndirectPhotonsWanted,
                                  indirectPhotons.size(), BlockSize))) {
                    Error("Unable to store enough photons.  Giving up.\n");
                    causticPhotons.erase(causticPhotons.begin(), causticPhotons.end());
                    indirectPhotons.erase(indirectPhotons.begin(), indirectPhotons.end());
                    radiancePhotons.erase(radiancePhotons.begin(), radiancePhotons.end());
                    abortTasks = true;
                    return;
                }

                int64_t numProgressed = 0;
                nShot += BlockSize;

                // Merge indirect photons into shared array
                if (!indirectDone) {
                    numProgressed += localIndirectPhotons.size();
                    this->nIndirectPaths += BlockSize;
                    for (uint32_t i = 0; i < localIndirectPhotons.size(); ++i)
                        indirectPhotons.push_back(localIndirectPhotons[i]);
                    localIndirectPhotons.erase(localIndirectPhotons.begin(),
                                               localIndirectPhotons.end());
                    if (indirectPhotons.size() >= this->nIndirectPhotonsWanted)
                        indirectDone = true;
                    nDirectPaths += BlockSize;
                    for (uint32_t i = 0; i < localDirectPhotons.size(); ++i)
                        directPhotons.push_back(localDirectPhotons[i]);
                    localDirectPhotons.erase(localDirectPhotons.begin(),
                                             localDirectPhotons.end());
                }

                // Merge direct, caustic, and radiance photons into shared array
                if (!causticDone) {
                    numProgressed += localCausticPhotons.size();
                    this->nCausticPaths += BlockSize;
                    for (uint32_t i = 0; i < localCausticPhotons.size(); ++i)
                        causticPhotons.push_back(localCausticPhotons[i]);
                    localCausticPhotons.erase(localCausticPhotons.begin(), localCausticPhotons.end());
                    if (causticPhotons.size() >= this->nCausticPhotonsWanted)
                        causticDone = true;
                }

                // Merge volume photons into shared array
                if (!volumeDone) {
                    numProgressed += localVolumePhotons.size();
                    for (uint32_t i = 0; i < localVolumePhotons.size(); ++i)
                        volumePhotons.push_back(localVolumePhotons[i]);
                    localVolumePhotons.erase(localVolumePhotons.begin(), localVolumePhotons.end());
                    if (volumePhotons.size() >= nVolumePhotonsWanted)
                        volumeDone = true;
                }

                progress.Update(numProgressed);

                for (uint32_t i = 0; i < localRadiancePhotons.size(); ++i)
                    radiancePhotons.push_back(localRadiancePhotons[i]);
                localRadiancePhotons.erase(localRadiancePhotons.begin(), localRadiancePhotons.end());
                for (uint32_t i = 0; i < localRpReflectances.size(); ++i)
                    rpReflectances.push_back(localRpReflectances[i]);
                localRpReflectances.erase(localRpReflectances.begin(), localRpReflectances.end());
                for (uint32_t i = 0; i < localRpTransmittances.size(); ++i)
                    rpTransmittances.push_back(localRpTransmittances[i]);
                localRpTransmittances.erase(localRpTransmittances.begin(), localRpTransmittances.end());
            }

            // Exit task if enough photons have been found
            if (indirectDone && causticDone && volumeDone)
                break;
        }

    }, nTasks);
    progress.Done();

    indirectPhotonCounter = indirectPhotons.size();
    causticPhotonCounter = causticPhotons.size();
    radiancePhotonCounter = radiancePhotons.size();
    volumePhotonCounter = volumePhotons.size();
}


void
PhotonIntegrator::BuildPhotonMaps(std::vector<Photon> const &directPhotons,
                                  std::vector<Photon> const &causticPhotons,
                                  std::vector<Photon> const &indirectPhotons,
                                  std::vector<Photon> const &volumePhotons) {
    ProgressReporter progress(3, "Building photon maps");
    ParallelFor([&](int64_t TaskIndex) {
        DCHECK_GE(TaskIndex, 0);
        DCHECK_LE(TaskIndex, 3);

        if (TaskIndex == 0 && directPhotons.size() > 0)
            directMap = new KdTree<Photon>(directPhotons);
        else if (TaskIndex == 1 && indirectPhotons.size() > 0)
            indirectMap = new KdTree<Photon>(indirectPhotons);
        else if (TaskIndex == 2 && causticPhotons.size() > 0)
            causticMap = new KdTree<Photon>(causticPhotons);
        else if (TaskIndex == 3 && volumePhotons.size() > 0)
            volumeMap = new KdTree<Photon>(volumePhotons);

        progress.Update(1);
    }, 4);
    progress.Done();
}


void PhotonIntegrator::ComputePhotonRadiances(std::vector<RadiancePhoton> &radiancePhotons,
                                              std::vector<Spectrum> const &rpReflectances,
                                              std::vector<Spectrum> const &rpTransmittances,
                                              KdTree<Photon> const *directMap,
                                              int const nDirectPaths) {
}


PhotonIntegrator *CreatePhotonMapIntegrator(
        const ParamSet &params, std::shared_ptr<Sampler> sampler,
        std::shared_ptr<const Camera> camera) {
    // TODO from pbrt-v2 PhotonIntegrator
    bool requirePhotons = params.FindOneBool("requirephotons", false);
    int nCaustic = params.FindOneInt("causticphotons", 20000);
    int nIndirect = params.FindOneInt("indirectphotons", 100000);
    int nVolume = params.FindOneInt("volumephotons", 50000);
    int nUsed = params.FindOneInt("nused", 50);
    if (PbrtOptions.quickRender) nCaustic = nCaustic / 10;
    if (PbrtOptions.quickRender) nIndirect = nIndirect / 10;
    if (PbrtOptions.quickRender) nUsed = std::max(1, nUsed / 10);
    int maxSpecularDepth = params.FindOneInt("maxspeculardepth", 5);
    int maxPhotonDepth = params.FindOneInt("maxphotondepth", 5);
    bool finalGather = params.FindOneBool("finalgather", true);
    int gatherSamples = params.FindOneInt("finalgathersamples", 32);
    if (PbrtOptions.quickRender) gatherSamples = std::max(1, gatherSamples / 4);
    Float maxDist = params.FindOneFloat("maxdist", .1f);
    Float gatherAngle = params.FindOneFloat("gatherangle", 10.f);

    // TODO maxDepth needed?
    // int maxDepth = params.FindOneInt("maxdepth", 5);
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]},
                                             {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }

    return new PhotonIntegrator(camera, sampler, pixelBounds,
                                nCaustic, nIndirect, nVolume, requirePhotons,
                                nUsed, maxSpecularDepth, maxPhotonDepth, maxDist, finalGather, gatherSamples,
                                gatherAngle);
}

template<size_t N>
std::array<Float, N> GetSamples1D(Sampler &sampler) {
    std::array<Float, N> values;
    for (int i = 0; i < N; ++i)
        values[i] = sampler.Get1D();
    return values;
};

template<size_t N>
std::array<Point2f, N> GetSamples2D(Sampler &sampler) {
    std::array<Point2f, N> values;
    for (int i = 0; i < N; ++i)
        values[i] = sampler.Get2D();
    return values;
};
}
