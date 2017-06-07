
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
#include "base/mutex.h"
#include "samplers/halton.h"
#include "integrators/photonmap.h"
#include "scene.h"
#include "progressreporter.h"
#include "paramset.h"
#include "camera.h"

namespace pbrt {
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

static Spectrum LPhoton(KdTree<Photon> *map, int nPaths, int nLookup,
                        ClosePhoton *lookupBuf, BSDF *bsdf, RNG &rng, const SurfaceInteraction &isect,
                        const Vector3f &w, Float maxDistSquared);

static Spectrum EPhoton(KdTree<Photon> *map, int count, int nLookup,
                        ClosePhoton *lookupBuf, Float maxDist2, const Point3f &p, const Normal3f &n);

// PhotonIntegrator Local Definitions
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


Spectrum LPhoton(KdTree<Photon> *map, int nPaths, int nLookup,
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
            int sqrtSamples = 6;
            int nSamples = sqrtSamples * sqrtSamples;
            Point2f const *s1 = sampler.Get2DArray(nSamples);
            Point2f const *s2 = sampler.Get2DArray(nSamples);

            L += Lr * bsdf->rho(wo, nSamples, s1, BxDFType(BSDF_ALL_TYPES | BSDF_REFLECTION)) * InvPi +
                 Lt * bsdf->rho(wo, nSamples, s2, BxDFType(BSDF_ALL_TYPES | BSDF_TRANSMISSION)) * InvPi;
        }
    }
    return L;
}


Spectrum EPhoton(KdTree<Photon> *map, int count, int nLookup,
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
    return E / (count * md2 * M_PI);
}


// PhotonIntegrator Method Definitions
PhotonIntegrator::PhotonIntegrator(std::shared_ptr<const Camera> camera,
                                   std::shared_ptr<Sampler> sampler,
                                   const Bounds2i &pixelBounds,
                                   int ncaus, int nind,
                                   int nl, int mdepth, int mphodepth, float mdist, bool fg,
                                   int gs, float ga) :
        SamplerIntegrator(camera, sampler, pixelBounds),
        nCausticPhotonsWanted(ncaus),
        nIndirectPhotonsWanted(nind),
        nLookup(nl),
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
    std::vector<Float> lightPbs(numLights);
    Float sumLightPowers = 0.0;
    for (int i = 0; i < numLights; ++i) {
        lightPbs[i] = scene.lights[i]->Power().MaxComponentValue();
        sumLightPowers += lightPbs[i];
    }
    DCHECK_GT(sumLightPowers, 0.0);
    for (int i = 0; i < numLights; ++i)
        lightPbs[i] /= sumLightPowers;
    lightDistr = std::unique_ptr<Distribution1D>(new Distribution1D(lightPbs.data(), numLights));

    if (scene.lights.size() == 0) return;

    // Declare shared variables for photon shooting
    Mutex mutex;
    int nDirectPaths = 0;
    std::vector<Photon> causticPhotons, directPhotons, indirectPhotons;
    std::vector<RadiancePhoton> radiancePhotons;
    bool abortTasks = false;
    causticPhotons.reserve(nCausticPhotonsWanted);
    indirectPhotons.reserve(nIndirectPhotonsWanted);
    uint32_t nshot = 0;
    std::vector<Spectrum> rpReflectances, rpTransmittances;

    // Run parallel tasks for photon shooting
    ProgressReporter progress(nCausticPhotonsWanted + nIndirectPhotonsWanted, "Shooting photons");
    int nTasks = NumSystemCores();
    ParallelFor([&](int i_task) {
        // shoot photons on this process
        // Loop until desired number of photons met:
        // 1. shoot n photons
        // 2. append local direct, indirect and caustic photon vectors to global vectors

        // Declare local variables for _PhotonShootingTask_
        MemoryArena arena;
        std::vector<Photon> localDirectPhotons, localIndirectPhotons, localCausticPhotons;
        std::vector<RadiancePhoton> localRadiancePhotons;
        uint32_t totalPaths = 0;
        bool causticDone = (this->nCausticPhotonsWanted == 0);
        bool indirectDone = (this->nIndirectPhotonsWanted == 0);
        std::vector<Spectrum> localRpReflectances, localRpTransmittances;

        while (true) {
            // Follow photon paths for a block of samples
            const uint32_t blockSize = 4096;
            for (uint32_t i = 0; i < blockSize; ++i) {
                Float const* u = sampler.Get1DArray(3);

                // Choose light to shoot photon from
                Float lightPdf;
                int lightNum = lightDistr->SampleDiscrete(u[0], &lightPdf);
                std::shared_ptr<Light> light = scene.lights[lightNum];

                // Generate _photonRay_ from light source and initialize _alpha_
                RayDifferential photonRay;
                Vector3f wi;
                Float pdf;
                VisibilityTester vis;
                Normal3f Nl;
                Interaction ref;
                Spectrum Le = light->Sample_Li(ref, Point2f(u[1], u[2]), &wi, &pdf, &vis);

                if (pdf == 0.f || Le.IsBlack()) continue;
                Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf);
                if (!alpha.IsBlack()) {
                    // Follow photon path through scene and record intersections
                    bool specularPath = true;
                    SurfaceInteraction photonIsect;
                    int nIntersections = 0;
                    while (scene.Intersect(photonRay, &photonIsect)) {
                        ++nIntersections;

                        // Handle photon/surface intersection
                        // TODO
                        // TODO
                        // TODO
                        // TODO
                        // TODO update alpha with volume integrator transmittance
                        // TODO => ray.medium.Tr(ray, sampler)
                        // TODO use isect.SpawnRay(direction) !!!!!!!!
                        // TODO
                        // TODO
                        // TODO
                        // TODO
                        // TODO
                        // TODO
                        // TODO
                        // TODO
                        // TODO
                        Spectrum transmittance(1.0f);
                        alpha *= transmittance;
                        BSDF *photonBSDF = photonIsect.bsdf;
                        BxDFType specularType = BxDFType(BSDF_REFLECTION |
                                                         BSDF_TRANSMISSION | BSDF_SPECULAR);
                        bool hasNonSpecular = (photonBSDF->NumComponents() >
                                               photonBSDF->NumComponents(specularType));
                        Vector3f wo = -photonRay.d;
                        if (hasNonSpecular) {
                            // Deposit photon at surface
                            Photon photon(photonIsect.p, alpha, wo);
                            bool depositedPhoton = false;
                            if (specularPath && nIntersections > 1) {
                                if (!causticDone) {
                                    depositedPhoton = true;
                                    localCausticPhotons.push_back(photon);
                                }
                            } else {
                                // Deposit either direct or indirect photon
                                // stop depositing direct photons once indirectDone is true; don't
                                // want to waste memory storing too many if we're going a long time
                                // trying to get enough caustic photons desposited.
                                if (nIntersections == 1 && !indirectDone && this->finalGather) {
                                    depositedPhoton = true;
                                    localDirectPhotons.push_back(photon);
                                } else if (nIntersections > 1 && !indirectDone) {
                                    depositedPhoton = true;
                                    localIndirectPhotons.push_back(photon);
                                }
                            }

                            // Possibly create radiance photon at photon intersection Point3f
                            if (depositedPhoton && this->finalGather &&
                                sampler.Get1D() < .125f) {
                                Normal3f n = photonIsect.n;
                                n = Faceforward(n, -photonRay.d);
                                localRadiancePhotons.push_back(RadiancePhoton(photonIsect.p, n));

                                int sqrtSamples = 6;
                                int nSamples = sqrtSamples * sqrtSamples;
                                Point2f const *s1 = sampler.Get2DArray(nSamples);
                                Point2f const *s2 = sampler.Get2DArray(nSamples);
                                Spectrum rho_r = photonBSDF->rho(nSamples, s1, s2, BxDFType(BSDF_ALL_TYPES | BSDF_REFLECTION));
                                localRpReflectances.push_back(rho_r);

                                s1 = sampler.Get2DArray(nSamples);
                                s2 = sampler.Get2DArray(nSamples);
                                Spectrum rho_t = photonBSDF->rho(nSamples, s1, s2, BxDFType(BSDF_ALL_TYPES | BSDF_TRANSMISSION));
                                localRpTransmittances.push_back(rho_t);
                            }
                        }
                        if (nIntersections >= this->maxPhotonDepth) break;

                        // Sample new photon ray direction
                        Vector3f wi;
                        Float pdf;
                        BxDFType flags;
                        Point2f uphoton = sampler.Get2D();
                        Spectrum fr = photonBSDF->Sample_f(wo, &wi, uphoton, &pdf, BSDF_ALL, &flags);
                        if (fr.IsBlack() || pdf == 0.f) break;
                        Spectrum anew = alpha * fr *
                                        AbsDot(wi, photonIsect.n) / pdf;

                        // Possibly terminate photon path with Russian roulette
                        Float continueProb = std::min(1.f, anew.y() / alpha.y());
                        if (sampler.Get1D() > continueProb)
                            break;
                        alpha = anew / continueProb;
                        specularPath &= ((flags & BSDF_SPECULAR) != 0);

                        if (indirectDone && !specularPath) break;
                        photonRay = RayDifferential(photonIsect.p, wi);
                    }
                }
                arena.Reset();
            }

            // Merge local photon data with data in _PhotonIntegrator_
            {
                MutexLock lock(&mutex);

                // Give up if we're not storing enough photons
                if (abortTasks)
                    return;
                if (nshot > 500000 &&
                    (unsuccessful(this->nCausticPhotonsWanted,
                                  causticPhotons.size(), blockSize) ||
                     unsuccessful(this->nIndirectPhotonsWanted,
                                  indirectPhotons.size(), blockSize))) {
                    Error("Unable to store enough photons.  Giving up.\n");
                    causticPhotons.erase(causticPhotons.begin(), causticPhotons.end());
                    indirectPhotons.erase(indirectPhotons.begin(), indirectPhotons.end());
                    radiancePhotons.erase(radiancePhotons.begin(), radiancePhotons.end());
                    abortTasks = true;
                    return;
                }
                progress.Update(localIndirectPhotons.size() + localCausticPhotons.size());
                nshot += blockSize;

                // Merge indirect photons into shared array
                if (!indirectDone) {
                    this->nIndirectPaths += blockSize;
                    for (uint32_t i = 0; i < localIndirectPhotons.size(); ++i)
                        indirectPhotons.push_back(localIndirectPhotons[i]);
                    localIndirectPhotons.erase(localIndirectPhotons.begin(),
                                               localIndirectPhotons.end());
                    if (indirectPhotons.size() >= this->nIndirectPhotonsWanted)
                        indirectDone = true;
                    nDirectPaths += blockSize;
                    for (uint32_t i = 0; i < localDirectPhotons.size(); ++i)
                        directPhotons.push_back(localDirectPhotons[i]);
                    localDirectPhotons.erase(localDirectPhotons.begin(),
                                             localDirectPhotons.end());
                }

                // Merge direct, caustic, and radiance photons into shared array
                if (!causticDone) {
                    this->nCausticPaths += blockSize;
                    for (uint32_t i = 0; i < localCausticPhotons.size(); ++i)
                        causticPhotons.push_back(localCausticPhotons[i]);
                    localCausticPhotons.erase(localCausticPhotons.begin(), localCausticPhotons.end());
                    if (causticPhotons.size() >= this->nCausticPhotonsWanted)
                        causticDone = true;
                }

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
            if (indirectDone && causticDone)
                break;
        }

    }, nTasks);
    progress.Done();

    // Build kd-trees for indirect and caustic photons
    KdTree<Photon> *directMap = nullptr;
    if (directPhotons.size() > 0)
        directMap = new KdTree<Photon>(directPhotons);
    if (causticPhotons.size() > 0)
        causticMap = new KdTree<Photon>(causticPhotons);
    if (indirectPhotons.size() > 0)
        indirectMap = new KdTree<Photon>(indirectPhotons);

    // Precompute radiance at a subset of the photons
    if (finalGather && radiancePhotons.size()) {
        // Launch tasks to compute photon radiances
        ProgressReporter progRadiance(nTasks, "Computing photon radiances");
        ParallelFor([&](uint32_t i_task) {
            // Compute range of radiance photons to process in task
            uint32_t taskSize = radiancePhotons.size() / nTasks;
            uint32_t excess = radiancePhotons.size() % nTasks;
            uint32_t rpStart = std::min(i_task, excess) * (taskSize + 1) +
                               std::max(0, (int) i_task - (int) excess) * taskSize;
            uint32_t rpEnd = rpStart + taskSize + (i_task < excess ? 1 : 0);
            if (i_task == nTasks - 1) DCHECK_EQ(rpEnd, radiancePhotons.size());
            ClosePhoton *lookupBuf = new ClosePhoton[nLookup];
            for (uint32_t i = rpStart; i < rpEnd; ++i) {
                // Compute radiance for radiance photon _i_
                RadiancePhoton &rp = radiancePhotons[i];
                const Spectrum &rho_r = rpReflectances[i], &rho_t = rpTransmittances[i];
                if (!rho_r.IsBlack()) {
                    // Accumulate outgoing radiance due to reflected irradiance
                    Spectrum E = EPhoton(directMap, nDirectPaths, nLookup, lookupBuf,
                                         maxDistSquared, rp.p, rp.n) +
                                 EPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                                         maxDistSquared, rp.p, rp.n) +
                                 EPhoton(causticMap, nCausticPaths, nLookup, lookupBuf,
                                         maxDistSquared, rp.p, rp.n);
                    rp.Lo += InvPi * rho_r * E;
                }
                if (!rho_t.IsBlack()) {
                    // Accumulate outgoing radiance due to transmitted irradiance
                    Spectrum E = EPhoton(directMap, nDirectPaths, nLookup, lookupBuf,
                                         maxDistSquared, rp.p, -rp.n) +
                                 EPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                                         maxDistSquared, rp.p, -rp.n) +
                                 EPhoton(causticMap, nCausticPaths, nLookup, lookupBuf,
                                         maxDistSquared, rp.p, -rp.n);
                    rp.Lo += InvPi * rho_t * E;
                }
            }
            delete[] lookupBuf;
            progress.Update();
        }, 64, 1);
        progRadiance.Done();
        radianceMap = new KdTree<RadiancePhoton>(radiancePhotons);
    }
    delete directMap;
}


Spectrum PhotonIntegrator::Li(const RayDifferential &ray,
                              const Scene &scene,
                              Sampler &sampler,
                              MemoryArena &arena,
                              int depth) const {
    ProfilePhase prof(Prof::SamplerIntegratorLi);

    Spectrum L(0.0f);

    // Find closest ray intersection or return background radiance
    SurfaceInteraction isect;
    if (!scene.Intersect(ray, &isect)) {
        for (const auto &light : scene.lights) L += light->Le(ray);
        return L;
    }

    Vector3f wo = -ray.d;

    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit Point3f
    BSDF *bsdf = isect.bsdf;
    const Point3f &p = isect.p;
    const Normal3f &n = isect.n;

    // TODO populate nLightSamples
    // TODO populate nLightSamples
    // TODO populate nLightSamples
    // TODO populate nLightSamples
    // TODO populate nLightSamples
    // TODO populate nLightSamples
    // TODO do what directLightingIntegrator does
    std::vector<int> nLightSamples;
    bool handleMedia = false;
    L += UniformSampleAllLights(isect, scene, arena, sampler, nLightSamples, handleMedia);
    // Compute caustic lighting for photon map integrator
    ClosePhoton *lookupBuf = arena.Alloc<ClosePhoton>(nLookup);
    L += LPhoton(causticMap, nCausticPaths, nLookup, lookupBuf, bsdf,
                 sampler, isect, wo, maxDistSquared);

    // Compute indirect lighting for photon map integrator
    if (finalGather && indirectMap != NULL) {
#if 1
        // Do one-bounce final gather for photon map
        BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
                                        BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
        if (bsdf->NumComponents(nonSpecular) > 0) {
            // Find indirect photons around Point3f  for importance sampling
            const uint32_t nIndirSamplePhotons = 50;
            PhotonProcess proc(nIndirSamplePhotons,
                               arena.Alloc<ClosePhoton>(nIndirSamplePhotons));
            Float searchDist2 = maxDistSquared;
            while (proc.nFound < nIndirSamplePhotons) {
                Float md2 = searchDist2;
                proc.nFound = 0;
                indirectMap->Lookup(p, proc, md2);
                searchDist2 *= 2.f;
            }

            // Copy photon directions to local array
            Vector3f *photonDirs = arena.Alloc<Vector3f>(nIndirSamplePhotons);
            for (uint32_t i = 0; i < nIndirSamplePhotons; ++i)
                photonDirs[i] = proc.photons[i].photon->wi;

            // Use BSDF to do final gathering
            Spectrum Li = 0.;
            for (int i = 0; i < gatherSamples; ++i) {
                // Sample random direction from BSDF for final gather ray
                Vector3f wi;
                Float pdf;
                Point2f u_bsdf;
                Spectrum fr = bsdf->Sample_f(wo, &wi, u_bsdf,
                                             &pdf, BxDFType(BSDF_ALL & ~BSDF_SPECULAR));
                if (fr.IsBlack() || pdf == 0.f) continue;
                DCHECK(pdf >= 0.f);

                // Trace BSDF final gather ray and accumulate radiance
                RayDifferential bounceRay(p, wi);
                SurfaceInteraction gatherIsect;
                if (scene.Intersect(bounceRay, &gatherIsect)) {
                    // Compute exitant radiance _Lindir_ using radiance photons
                    Spectrum Lindir = 0.f;
                    Normal3f nGather = gatherIsect.n;
                    nGather = Faceforward(nGather, -bounceRay.d);
                    RadiancePhotonProcess proc(nGather);
                    Float md2 = INFINITY;
                    radianceMap->Lookup(gatherIsect.p, proc, md2);
                    if (proc.photon != NULL)
                        Lindir = proc.photon->Lo;

                    // TODO compute transmittance
                    // TODO compute transmittance
                    // TODO compute transmittance
                    // TODO compute transmittance
                    // TODO compute transmittance
                    // TODO compute transmittance
                    Spectrum transmittance(1.0f);
                    Lindir *= transmittance;

                    // Compute MIS weight for BSDF-sampled gather ray

                    // Compute PDF for photon-sampling of direction _wi_
                    Float photonPdf = 0.f;
                    Float conePdf = UniformConePdf(cosGatherAngle);
                    for (uint32_t j = 0; j < nIndirSamplePhotons; ++j)
                        if (Dot(photonDirs[j], wi) > .999f * cosGatherAngle)
                            photonPdf += conePdf;
                    photonPdf /= nIndirSamplePhotons;
                    Float wt = PowerHeuristic(gatherSamples, pdf, gatherSamples, photonPdf);
                    Li += fr * Lindir * (AbsDot(wi, n) * wt / pdf);
                }
            }
            L += Li / gatherSamples;

            // Use nearby photons to do final gathering
            Li = 0.;
            for (int i = 0; i < gatherSamples; ++i) {
                // Sample random direction using photons for final gather ray
                // TODO
                // TODO
                // TODO
                // TODO
                // TODO what is BSDFSample equivalent?
                // TODO => random sampler.2D
                // TODO
                // TODO
                // TODO
                Point2f u_cone = sampler.Get2D();
                int photonNum = std::min((int) nIndirSamplePhotons - 1,
                                         static_cast<int>(std::floor(sampler.Get1D() * nIndirSamplePhotons)));

                // Sample gather ray direction from _photonNum_
                Vector3f vx, vy;
                CoordinateSystem(photonDirs[photonNum], &vx, &vy);
                Vector3f wi = UniformSampleCone(u_cone, cosGatherAngle, vx, vy, photonDirs[photonNum]);
                // Trace photon-sampled final gather ray and accumulate radiance
                Spectrum fr = bsdf->f(wo, wi);
                if (fr.IsBlack()) continue;
                RayDifferential bounceRay(p, wi);
                SurfaceInteraction gatherIsect;
                if (scene.Intersect(bounceRay, &gatherIsect)) {
                    // Compute exitant radiance _Lindir_ using radiance photons
                    Spectrum Lindir = 0.f;
                    Normal3f nGather = gatherIsect.n;
                    nGather = Faceforward(nGather, -bounceRay.d);
                    RadiancePhotonProcess proc(nGather);
                    Float md2 = INFINITY;
                    radianceMap->Lookup(gatherIsect.p, proc, md2);
                    if (proc.photon != NULL)
                        Lindir = proc.photon->Lo;

                    // TODO compute transmittance
                    // TODO compute transmittance
                    // TODO compute transmittance
                    // TODO compute transmittance
                    // TODO compute transmittance
                    // TODO compute transmittance
                    Spectrum transmittance(1.0f);
                    Lindir *= transmittance;

                    // Compute PDF for photon-sampling of direction _wi_
                    Float photonPdf = 0.f;
                    Float conePdf = UniformConePdf(cosGatherAngle);
                    for (uint32_t j = 0; j < nIndirSamplePhotons; ++j)
                        if (Dot(photonDirs[j], wi) > .999f * cosGatherAngle)
                            photonPdf += conePdf;
                    photonPdf /= nIndirSamplePhotons;

                    // Compute MIS weight for photon-sampled gather ray
                    Float bsdfPdf = bsdf->Pdf(wo, wi);
                    Float wt = PowerHeuristic(gatherSamples, photonPdf, gatherSamples, bsdfPdf);
                    Li += fr * Lindir * AbsDot(wi, n) * wt / photonPdf;
                }
            }
            L += Li / gatherSamples;
        }
#else
        // for debugging / examples: use the photon map directly
        Normal3f nn = Faceforward(n, -ray.d);
        RadiancePhotonProcess proc(nn);
        Float md2 = INFINITY;
        radianceMap->Lookup(p, proc, md2);
        if (proc.photon)
            L += proc.photon->Lo;
#endif
    } else
        L += LPhoton(indirectMap, nIndirectPaths, nLookup, lookupBuf,
                     bsdf, sampler, isect, wo, maxDistSquared);
    if (depth + 1 < maxSpecularDepth) {
        Vector3f wi;

        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, isect, scene, sampler, arena, depth);
        L += SpecularTransmit(ray, isect, scene, sampler, arena, depth);
    }
    return L;
}


PhotonIntegrator *CreatePhotonMapSurfaceIntegrator(
        const ParamSet &params, std::shared_ptr<Sampler> sampler,
        std::shared_ptr<const Camera> camera) {
    // TODO from pbrt-v2 PhotonIntegrator
    int nCaustic = params.FindOneInt("causticphotons", 20000);
    int nIndirect = params.FindOneInt("indirectphotons", 100000);
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
                                nCaustic, nIndirect,
                                nUsed, maxSpecularDepth, maxPhotonDepth, maxDist, finalGather, gatherSamples,
                                gatherAngle);
}
}
