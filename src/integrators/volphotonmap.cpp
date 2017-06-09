/**
 * Created by Benjamin Wiberg on 08/06/2017.
 */

#include "volphotonmap.h"
#include <array>
#include <samplers/zerotwosequence.h>
#include <base/mutex.h> // from glog
#include <samplers/halton.h>
#include "scene.h"
#include "progressreporter.h"
#include "paramset.h"
#include "camera.h"

namespace pbrt {
//STAT_COUNTER("Photon Mapping/Total Photon Surface Interactions ", totalPhotonSurfaceInteractions);
//STAT_COUNTER("Photon Mapping/Total Photon Volume Interactions ", totalPhotonVolumeInteractions);
//STAT_COUNTER("Photon Mapping/Photon paths followed", photonPaths);
//STAT_COUNTER("Photon Mapping/Total indirect photons", indirectPhotonCounter);
//STAT_COUNTER("Photon Mapping/Total caustic photons", causticPhotonCounter);
//STAT_COUNTER("Photon Mapping/Total radiance photons", radiancePhotonCounter);
//STAT_COUNTER("Photon Mapping/Total volume photons", volumePhotonCounter);
//STAT_MEMORY_COUNTER("Memory/Photon maps", photonMapBytes);
//
//struct Photon {
//    Photon(const Point3f &pp, const Spectrum &wt, const Vector3f &w)
//            : p(pp), alpha(wt), wi(w) {}
//
//    Photon() {}
//
//    Point3f p;
//    Spectrum alpha;
//    Vector3f wi;
//};
//
//struct RadiancePhoton {
//    RadiancePhoton(const Point3f &pp, const Normal3f &nn)
//            : p(pp), n(nn), Lo(0.f) {}
//
//    RadiancePhoton() {}
//
//    Point3f p;
//    Normal3f n;
//    Spectrum Lo;
//};
//
//
//struct PhotonProcess {
//    // PhotonProcess Public Methods
//    PhotonProcess(uint32_t mp, ClosePhoton *buf);
//
//    void operator()(const Point3f &p, const Photon &photon, Float dist2,
//                    Float &maxDistSquared);
//
//    ClosePhoton *photons;
//    uint32_t nLookup, nFound;
//};
//
//
//struct ClosePhoton {
//    // ClosePhoton Public Methods
//    ClosePhoton(const Photon *p = NULL, Float md2 = INFINITY)
//            : photon(p), distanceSquared(md2) {}
//
//    bool operator<(const ClosePhoton &p2) const {
//        return distanceSquared == p2.distanceSquared ?
//               (photon < p2.photon) : (distanceSquared < p2.distanceSquared);
//    }
//
//    const Photon *photon;
//    Float distanceSquared;
//};
//
//
//PhotonProcess::PhotonProcess(uint32_t mp, ClosePhoton *buf)
//        : photons(buf), nLookup(mp), nFound(0) {}
//
//
//struct RadiancePhotonProcess {
//    // RadiancePhotonProcess Methods
//    RadiancePhotonProcess(const Normal3f &nn)
//            : n(nn) {
//        photon = NULL;
//    }
//
//    void operator()(const Point3f &p, const RadiancePhoton &rp,
//                    Float distSquared, Float &maxDistSquared) {
//        if (Dot(rp.n, n) > 0) {
//            photon = &rp;
//            maxDistSquared = distSquared;
//        }
//    }
//
//    const Normal3f &n;
//    const RadiancePhoton *photon;
//};
//
//
//inline Float kernel(const Photon *photon, const Point3f &p, Float maxDist2);
//
//static Spectrum LPhoton(KdTree<Photon> const *map, int nPaths, int nLookup,
//                        ClosePhoton *lookupBuf, BSDF *bsdf, RNG &rng, const SurfaceInteraction &isect,
//                        const Vector3f &w, Float maxDistSquared);
//
//static Spectrum EPhoton(KdTree<Photon> const *map, int count, int nLookup,
//                        ClosePhoton *lookupBuf, Float maxDist2, const Point3f &p, const Normal3f &n);
//
//// VolPhotonIntegrator Local Deff
//// ifinitions
//inline bool unsuccessful(uint32_t needed, uint32_t found, uint32_t shot) {
//    return (found < needed && (found == 0 || found < shot / 1024));
//}
//
//
//inline void PhotonProcess::operator()(const Point3f &p,
//                                      const Photon &photon, Float distSquared, Float &maxDistSquared) {
//    if (nFound < nLookup) {
//        // Add photon to unordered array of photons
//        photons[nFound++] = ClosePhoton(&photon, distSquared);
//        if (nFound == nLookup) {
//            std::make_heap(&photons[0], &photons[nLookup]);
//            maxDistSquared = photons[0].distanceSquared;
//        }
//    } else {
//        // Remove most distant photon from heap and add new photon
//        std::pop_heap(&photons[0], &photons[nLookup]);
//        photons[nLookup - 1] = ClosePhoton(&photon, distSquared);
//        std::push_heap(&photons[0], &photons[nLookup]);
//        maxDistSquared = photons[0].distanceSquared;
//    }
//}
//
//
//inline Float kernel(const Photon *photon, const Point3f &p,
//                    Float maxDist2) {
//    Float s = (1.f - DistanceSquared(photon->p, p) / maxDist2);
//    return 3.f * InvPi * s * s;
//}
//
//
//Spectrum LPhoton(KdTree<Photon> const *map, int nPaths, int nLookup,
//                 ClosePhoton *lookupBuf, BSDF *bsdf, Sampler &sampler,
//                 const SurfaceInteraction &isect, const Vector3f &wo, Float maxDist2) {
//    return Spectrum(0.f);
//}
//
//
//Spectrum EPhoton(KdTree<Photon> const *map, int count, int nLookup,
//                 ClosePhoton *lookupBuf, Float maxDist2, const Point3f &p,
//                 const Normal3f &n) {
//    return Spectrum(0.f);
//}


// VolPhotonIntegrator Method Definitions
VolPhotonIntegrator::VolPhotonIntegrator(std::shared_ptr<const Camera> camera,
                                         std::shared_ptr<Sampler> sampler,
                                         const Bounds2i &pixelBounds,
                                         int ncaus, int nind, int nvol, bool reqphotons,
                                         int nl,
                                         int mspecdepth, int mphodepth, int samplesperdepth,
                                         float mdist, bool fg,
                                         int gs, float ga) :
        SamplerIntegrator(camera, sampler, pixelBounds),
        nCausticPhotonsWanted(ncaus),
        nIndirectPhotonsWanted(nind),
        nVolumePhotonsWanted(nvol),
        nLookup(nl),
        requirePhotons(reqphotons),
        maxSpecularDepth(mspecdepth),
        maxPhotonDepth(mphodepth),
        samplesPerPhotonDepth(samplesperdepth),
        maxDistSquared(mdist * mdist),
        finalGather(fg),
        cosGatherAngle(cos(Radians(ga))),
        gatherSamples(gs),
        nCausticPaths(nIndirectPaths = 0),
        volumeMap(nullptr) {}


VolPhotonIntegrator::~VolPhotonIntegrator() {
    delete volumeMap;
}


void VolPhotonIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
//    if (scene.lights.size() == 0) return;
//
//    ComputeLightSamples(scene, sampler);
//    ComputeLightPhotonDistrib(scene);
//
//    // Declare shared variables for photon shooting processes
//    std::vector<Photon> directPhotons, indirectPhotons, causticPhotons, volumePhotons;
//    directPhotons.reserve(64 * nIndirectPhotonsWanted);
//    indirectPhotons.reserve(nIndirectPhotonsWanted);
//    causticPhotons.reserve(nCausticPhotonsWanted);
//    volumePhotons.reserve(nVolumePhotonsWanted);
//    std::vector<RadiancePhoton> radiancePhotons;
//    std::vector<Spectrum> radpReflectances, radpTransmittances;
//    uint32_t nDirectPaths = 0, nShot = 0;
//
//    ShootPhotons(scene,
//                 directPhotons, indirectPhotons, causticPhotons, volumePhotons,
//                 radiancePhotons, radpReflectances, radpTransmittances,
//                 nDirectPaths, nShot);
//    BuildPhotonMaps(directPhotons, indirectPhotons, causticPhotons, volumePhotons);
}


Spectrum VolPhotonIntegrator::Li(const RayDifferential &ray,
                                 const Scene &scene,
                                 Sampler &sampler,
                                 MemoryArena &arena,
                                 int depth) const {
    Spectrum L(0.0f);

//    // Find closest ray intersection or return background radiance
//    SurfaceInteraction isect;
//    if (!scene.Intersect(ray, &isect)) {
//        for (const auto &light : scene.lights) L += light->Le(ray);
//        return L;
//    }
//
//    const Point3f p = isect.p;
//
//    const uint32_t nIndirSamplePhotons = 50;
//    PhotonProcess proc(nIndirSamplePhotons,
//                       arena.Alloc<ClosePhoton>(nIndirSamplePhotons));
//    Float md2 = maxDistSquared;
//    indirectMap->Lookup(p, proc, md2);
//    if (proc.nFound > 0)
//        L = Spectrum(1.0f);

    return L;
}


void VolPhotonIntegrator::ComputeLightSamples(const Scene &scene, Sampler &sampler) {
    for (int i = 0; i < scene.lights.size(); ++i) {
        // Compute number of samples to use for each light
        nLightSamples.push_back(sampler.RoundCount(scene.lights[i]->nSamples));
    }
}


void VolPhotonIntegrator::ComputeLightPhotonDistrib(Scene const &scene) {
    // compute power-based probability distribution of lights
    auto const numLights = scene.lights.size();

    std::vector<Float> lightPbs(numLights);
    Float sumLightPowers = 0.0;
    for (int i = 0; i < numLights; ++i) {
        lightPbs[i] = scene.lights[i]->Power().y();
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
VolPhotonIntegrator::ShootPhotons(Scene const &scene,
                                  std::vector<Photon> &directPhotons,
                                  std::vector<Photon> &indirectPhotons,
                                  std::vector<Photon> &causticPhotons,
                                  std::vector<Photon> &volumePhotons,
                                  std::vector<RadiancePhoton> &radiancePhotons,
                                  std::vector<Spectrum> &rpReflectances, std::vector<Spectrum> &rpTransmittances,
                                  uint32_t &nDirectPaths, uint32_t &nShot) {
//    ProgressReporter progress(nCausticPhotonsWanted + nCausticPhotonsWanted + nVolumePhotonsWanted,
//                              "Shooting photons");
//    bool abortTasks = false;
//
//    const int NumTasks = 1; //NumSystemCores();
//    Mutex mutex;
//
//    uint32_t totalPhotonsSent = 0;
//
//    ParallelFor([&](const int64_t TaskIndex) {
//        MemoryArena arena;
//
//        constexpr uint32_t BatchSize = 512;
//        std::vector<Photon> localDirectPhotons, localIndirectPhotons, localCausticPhotons, localVolumePhotons;
//        localDirectPhotons.reserve(BatchSize);
//        localIndirectPhotons.reserve(BatchSize);
//        localCausticPhotons.reserve(BatchSize);
//        localVolumePhotons.reserve(BatchSize);
//        std::vector<RadiancePhoton> localRadiancePhotons;
//        std::vector<Spectrum> localRadpRfls, localRadpTrs;
//        uint32_t totalPaths = 0;
//
//        bool indirectDone = nIndirectPhotonsWanted == 0;
//        bool causticDone = nCausticPhotonsWanted == 0;
//        bool volumeDone = nVolumePhotonsWanted == 0;
//
//        // Shoot photons until required number of photons stored
//        while (true) {
//            uint32_t LocalPhotonStartIndex;
//            {
//                MutexLock lock(&mutex);
//                LocalPhotonStartIndex = totalPhotonsSent;
//                totalPhotonsSent += BatchSize;
//            }
//
//            for (uint32_t i = 0; i < BatchSize; ++i) {
//                const uint32_t PhotonIndex = LocalPhotonStartIndex + i;
//                RNG rng(PhotonIndex);
//                AwesomeSampler sampler(0, PhotonIndex, rng);
//
//                auto Get2Ds = [&](const size_t NumValues) {
//                    Point2f *values = arena.Alloc<Point2f>(NumValues, false);
//                    for (size_t j = 0; j < NumValues; ++j)
//                        values[j] = sampler.Get2D();
//                    return values;
//                };
//
//
//                // Choose light to shoot photon from
//                Float lightPdf;
//                const int LightNum = lightDistr->SampleDiscrete(sampler.Get1D(), &lightPdf);
//                std::shared_ptr<Light> light = scene.lights[LightNum];
//
//                // Generate photon ray from light source and init alph
//                RayDifferential photonRay;
//                Float pdfPos, pdfDir;
//                Normal3f lightNormal;
//                Spectrum Le = light->Sample_Le(sampler.Get2D(), sampler.Get2D(),
//                                               0.0f, // hardcode time to zero, no moving things
//                                               &photonRay, &lightNormal, &pdfPos, &pdfDir);
//                if (pdfPos == 0.0f || pdfDir == 0.0f || Le.IsBlack()) continue;
//                ++photonPaths;
//                Spectrum alpha = (AbsDot(lightNormal, photonRay.d) * Le) / (pdfPos * pdfDir * lightPdf);
//                if (alpha.IsBlack()) continue;
//
//                // WIIE
//                // Follow photon through scene and record intersections
//                bool isSpecPath = true;
//                SurfaceInteraction photonIsect;
//                uint32_t nIntersections = 0;
//                while (scene.Intersect(photonRay, &photonIsect)) {
//                    ++nIntersections;
//
//                    // Participating media interactions
//                    MediumInteraction mi;
//                    if (photonRay.medium) alpha *= photonRay.medium->Sample(photonRay, sampler, arena, &mi);
//                    if (alpha.IsBlack()) break;
//
//                    // Handle photon interaction with medium OR surface
//                    photonIsect.ComputeScatteringFunctions(photonRay, arena);
//                    BSDF *bsdf = photonIsect.bsdf;
//                    if (!bsdf) {
//                        photonRay = photonIsect.SpawnRay(photonRay.d);
//                        nIntersections--;
//                        continue;
//                    }
//
//                    Vector3f wo = -photonRay.d;
//                    if (mi.IsValid()) {
//                        ++totalPhotonVolumeInteractions;
//
//                        Photon photon(mi.p, alpha, wo);
//                        if (nIntersections > 1 && !volumeDone)
//                            localVolumePhotons.push_back(photon);
//
//                    } else {
//                        ++totalPhotonSurfaceInteractions;
//
//                        // Handle photon/surface interaction
//                        constexpr BxDFType SpecularBxDFType = BxDFType(BSDF_REFLECTION |
//                                                                       BSDF_TRANSMISSION |
//                                                                       BSDF_SPECULAR);
//                        const bool HasNonSpecular = bsdf->NumComponents() >
//                                                    bsdf->NumComponents(SpecularBxDFType);
//                        if (HasNonSpecular) {
//                            // Deposit photon at surface
//                            const Photon photon(photonIsect.p, alpha, wo);
//                            bool hasDepositedPhoton = false;
//                            if (isSpecPath && nIntersections > 1 && !causticDone) {
//                                hasDepositedPhoton = true;
//                                localCausticPhotons.push_back(photon);
//                            } else {
//                                // Deposit either direct or indirect photon
//                                // stop depositing direct photons once indirectDone is true; don't
//                                // want to waste memory storing too many if we're going a long time
//                                // trying to get enough caustic photons desposited.
//                                if (nIntersections == 1 && !indirectDone && finalGather) {
//                                    hasDepositedPhoton = true;
//                                    localDirectPhotons.push_back(photon);
//                                } else if (nIntersections > 1 && !indirectDone) {
//                                    hasDepositedPhoton = true;
//                                    localIndirectPhotons.push_back(photon);
//                                }
//                            }
//
//                            // Possibly create radiance photon at photon intersection point
////                            if (hasDepositedPhoton && finalGather && sampler.Get1D() < .125f) {
////                                const Normal3f normal = Faceforward(photonIsect.n, -photonRay.d);
////                                localRadiancePhotons.push_back(RadiancePhoton(photonIsect.p, normal));
////
////                                constexpr int32_t SqrtSamples = 6;
////                                constexpr int32_t NumSamples = SqrtSamples * SqrtSamples;
////
////                                Point2f *s1 = Get2Ds(NumSamples);
////                                Point2f *s2 = Get2Ds(NumSamples);
////                                //std::array<Point2f, NumSamples> s1 = GetSamples2D<NumSamples>(*sampler);
////                                //std::array<Point2f, NumSamples> s2 = GetSamples2D<NumSamples>(*sampler);
////                                constexpr BxDFType ReflectanceBxDFType = BxDFType(BSDF_ALL_TYPES | BSDF_REFLECTION);
////                                localRadpRfls.push_back(bsdf->rho(NumSamples, s1, s2, ReflectanceBxDFType));
////
////                                s1 = Get2Ds(NumSamples);
////                                s2 = Get2Ds(NumSamples);
////                                constexpr BxDFType TransmittanceBxDFType = BxDFType(BSDF_ALL_TYPES | BSDF_TRANSMISSION);
////                                localRadpTrs.push_back(bsdf->rho(NumSamples, s1, s2, TransmittanceBxDFType));
////                            }
//                        }
//                    }
//
//                    if (nIntersections >= maxPhotonDepth) break;
//
//                    // Sample new photon ray direction
//                    Vector3f wi;
//                    Float pdf;
//                    BxDFType sampledBxDFType;
//                    Spectrum fr = bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, BSDF_ALL, &sampledBxDFType);
//
//                    if (pdf == 0.0f || fr.IsBlack()) break;
//                    Spectrum anew = alpha * fr * AbsDot(wi, photonIsect.n) / pdf;
//
//                    // Possibly terminate photon path with Russian roulette
//                    Float continueProb = std::min(1.0f, anew.y() / alpha.y());
//                    if (sampler.Get1D() > continueProb) break;
//
//                    alpha = anew / continueProb;
//                    isSpecPath &= (sampledBxDFType & BSDF_SPECULAR) != 0;
//
//                    if (indirectDone && !isSpecPath) break;
//                    photonRay = photonIsect.SpawnRay(wi);
//
//                    arena.Reset();
//                }
//                arena.Reset();
//            }
//
//            // Merge local photon data with data in _PhotonIntegrator_
//            {
//                MutexLock lock(&mutex);
//
//                if (abortTasks)
//                    return;
//                // Give up if we're not storing enough photons
//                if (!requirePhotons && nShot > 500000 &&
//                    (unsuccessful(nCausticPhotonsWanted,
//                                  causticPhotons.size(), BatchSize) ||
//                     unsuccessful(nIndirectPhotonsWanted,
//                                  indirectPhotons.size(), BatchSize))) {
//                    Error("Unable to store enough photons.  Giving up.\n");
//                    causticPhotons.erase(causticPhotons.begin(), causticPhotons.end());
//                    indirectPhotons.erase(indirectPhotons.begin(), indirectPhotons.end());
//                    radiancePhotons.erase(radiancePhotons.begin(), radiancePhotons.end());
//                    abortTasks = true;
//                    return;
//                }
//
//                progress.Update(localIndirectPhotons.size() + localCausticPhotons.size() + localVolumePhotons.size());
//                nShot += BatchSize;
//
//                // Merge indirect photons into shared array
//                if (!indirectDone) {
//                    nIndirectPaths += BatchSize;
//                    for (uint32_t i = 0; i < localIndirectPhotons.size(); ++i)
//                        indirectPhotons.push_back(localIndirectPhotons[i]);
//                    localIndirectPhotons.erase(localIndirectPhotons.begin(),
//                                               localIndirectPhotons.end());
//                    if (indirectPhotons.size() >= nIndirectPhotonsWanted)
//                        indirectDone = true;
//                    nDirectPaths += BatchSize;
//                    for (uint32_t i = 0; i < localDirectPhotons.size(); ++i)
//                        directPhotons.push_back(localDirectPhotons[i]);
//                    localDirectPhotons.erase(localDirectPhotons.begin(),
//                                             localDirectPhotons.end());
//                }
//
//                // Merge direct, caustic, and radiance photons into shared array
//                if (!causticDone) {
//                    nCausticPaths += BatchSize;
//                    for (uint32_t i = 0; i < localCausticPhotons.size(); ++i)
//                        causticPhotons.push_back(localCausticPhotons[i]);
//                    localCausticPhotons.erase(localCausticPhotons.begin(), localCausticPhotons.end());
//                    if (causticPhotons.size() >= nCausticPhotonsWanted)
//                        causticDone = true;
//                }
//
//                // Merge volume photons into shared array
//                if (!volumeDone) {
//                    for (uint32_t i = 0; i < localVolumePhotons.size(); ++i)
//                        volumePhotons.push_back(localVolumePhotons[i]);
//                    localVolumePhotons.erase(localVolumePhotons.begin(), localVolumePhotons.end());
//                    if (volumePhotons.size() >= nVolumePhotonsWanted)
//                        volumeDone = true;
//                }
//
//                for (uint32_t i = 0; i < localRadiancePhotons.size(); ++i)
//                    radiancePhotons.push_back(localRadiancePhotons[i]);
//                localRadiancePhotons.erase(localRadiancePhotons.begin(), localRadiancePhotons.end());
//                for (uint32_t i = 0; i < localRadpRfls.size(); ++i)
//                    rpReflectances.push_back(localRadpRfls[i]);
//                localRadpRfls.erase(localRadpRfls.begin(), localRadpRfls.end());
//                for (uint32_t i = 0; i < localRadpTrs.size(); ++i)
//                    rpTransmittances.push_back(localRadpTrs[i]);
//                localRadpTrs.erase(localRadpTrs.begin(), localRadpTrs.end());
//            }
//
//            // Exit task if enough photons have been found
//            if (indirectDone && causticDone)
//                break;
//        }
//    }, NumTasks);
//    progress.Done();
//
//    if (requirePhotons) {
//        DCHECK_GE(indirectPhotons.size(), nIndirectPhotonsWanted);
//        DCHECK_GE(causticPhotons.size(), nCausticPhotonsWanted);
//        DCHECK_GE(volumePhotons.size(), nVolumePhotonsWanted);
//    }
//
//    indirectPhotonCounter = indirectPhotons.size();
//    causticPhotonCounter = causticPhotons.size();
//    radiancePhotonCounter = radiancePhotons.size();
//    volumePhotonCounter = volumePhotons.size();
}


void VolPhotonIntegrator::BuildPhotonMaps(std::vector<Photon> const &directPhotons,
                                          std::vector<Photon> const &indirectPhotons,
                                          std::vector<Photon> const &causticPhotons,
                                          std::vector<Photon> const &volumePhotons) {
//    ProgressReporter progress(3, "Building photon maps");
//    ParallelFor([&](int64_t TaskIndex) {
//        DCHECK_GE(TaskIndex, 0);
//        DCHECK_LE(TaskIndex, 3);
//
//        if (TaskIndex == 0 && directPhotons.size() > 0)
//            directMap = new KdTree<Photon>(directPhotons);
//        else if (TaskIndex == 1 && indirectPhotons.size() > 0)
//            indirectMap = new KdTree<Photon>(indirectPhotons);
//        else if (TaskIndex == 2 && causticPhotons.size() > 0)
//            causticMap = new KdTree<Photon>(causticPhotons);
//        else if (TaskIndex == 3 && volumePhotons.size() > 0)
//            volumeMap = new KdTree<Photon>(volumePhotons);
//
//        progress.Update(1);
//    }, 4);
//    progress.Done();
}


void VolPhotonIntegrator::ComputePhotonRadiances(std::vector<RadiancePhoton> &radiancePhotons,
                                                 std::vector<Spectrum> const &rpReflectances,
                                                 std::vector<Spectrum> const &rpTransmittances,
                                                 KdTree<Photon> const *directMap,
                                                 int const nDirectPaths) {

}


VolPhotonIntegrator *CreateVolPhotonMapIntegrator(
        const ParamSet &params, std::shared_ptr<Sampler> sampler,
        std::shared_ptr<const Camera> camera) {
    bool requirePhotons = params.FindOneBool("requirephotons", false);
    int nCaustic = params.FindOneInt("causticphotons", 20000);
    int nIndirect = params.FindOneInt("indirectphotons", 100000);
    int nVolume = params.FindOneInt("volumephotons", 0); // TODO set some value
    int nUsed = params.FindOneInt("nused", 50);
    if (PbrtOptions.quickRender) nCaustic = nCaustic / 10;
    if (PbrtOptions.quickRender) nIndirect = nIndirect / 10;
    if (PbrtOptions.quickRender) nUsed = std::max(1, nUsed / 10);
    int maxSpecularDepth = params.FindOneInt("maxspeculardepth", 5);
    int maxPhotonDepth = params.FindOneInt("maxphotondepth", 5);
    int samplesPerPhotonDepth = params.FindOneInt("samplesperdepth", 240);
    bool finalGather = params.FindOneBool("finalgather", true);
    int gatherSamples = params.FindOneInt("finalgathersamples", 32);
    if (PbrtOptions.quickRender) gatherSamples = std::max(1, gatherSamples / 4);
    Float maxDist = params.FindOneFloat("maxdist", .1f);
    Float gatherAngle = params.FindOneFloat("gatherangle", 10.f);

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

    return new VolPhotonIntegrator(camera, sampler, pixelBounds,
                                   nCaustic, nIndirect, nVolume, requirePhotons,
                                   nUsed,
                                   maxSpecularDepth, maxPhotonDepth, samplesPerPhotonDepth,
                                   maxDist, finalGather, gatherSamples,
                                   gatherAngle);
}
}
