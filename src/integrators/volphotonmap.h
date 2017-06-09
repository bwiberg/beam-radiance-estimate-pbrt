
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_VOLPHOTONMAP_H
#define PBRT_INTEGRATORS_VOLPHOTONMAP_H

// integrators/photonmap.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"
#include "lightdistrib.h"

namespace pbrt {
struct Photon;
struct RadiancePhoton;
struct ClosePhoton;
struct PhotonProcess;
struct RadiancePhotonProcess;

// VolPhotonIntegrator Declarations
class VolPhotonIntegrator : public SamplerIntegrator {
public:
    // VolPhotonIntegrator Public Methods
    VolPhotonIntegrator(std::shared_ptr<const Camera> camera,
                        std::shared_ptr<Sampler> sampler,
                        const Bounds2i &pixelBounds,
                        int ncaus, int nindir, int nvol, bool reqphotons,
                        int nLookup,
                        int maxspecdepth, int maxphotondepth, int samplerperdepth,
                        float maxdist, bool finalGather, int gatherSamples,
                        float ga);

    ~VolPhotonIntegrator();

    Spectrum Li(const RayDifferential &ray,
                const Scene &scene,
                Sampler &sampler,
                MemoryArena &arena,
                int depth = 0) const;



    void Preprocess(const Scene &scene, Sampler &sampler);

private:
    void ComputeLightSamples(Scene const &scene, Sampler &sampler);

    void ComputeLightPhotonDistrib(Scene const &scene);

    void ShootPhotons(Scene const &scene,
                      std::vector<Photon> &directPhotons,
                      std::vector<Photon> &indirectPhotons,
                      std::vector<Photon> &causticPhotons,
                      std::vector<Photon> &volumePhotons,
                      std::vector<RadiancePhoton> &radiancePhotons,
                      std::vector<Spectrum> &rpReflectances, std::vector<Spectrum> &rpTransmittances,
                      uint32_t &nDirectPaths, uint32_t &nShot);

    void BuildPhotonMaps(std::vector<Photon> const &directPhotons,
                         std::vector<Photon> const &indirectPhotons,
                         std::vector<Photon> const &causticPhotons,
                         std::vector<Photon> const &volumePhotons);

    // VolPhotonIntegrator Private Methods
    void ComputePhotonRadiances(std::vector<RadiancePhoton> &radiancePhotons,
                                std::vector<Spectrum> const &rpReflectances,
                                std::vector<Spectrum> const &rpTransmittances,
                                KdTree<Photon> const *directMap,
                                int nDirectPaths);

    // VolPhotonIntegrator Private Data
    uint32_t nCausticPhotonsWanted, nIndirectPhotonsWanted, nVolumePhotonsWanted, nLookup;
    bool requirePhotons;
    Float maxDistSquared;

    int maxSpecularDepth, maxPhotonDepth;
    int samplesPerPhotonDepth;

    bool finalGather;
    int gatherSamples;
    Float cosGatherAngle;

    const std::string lightSampleStrategy;
    std::unique_ptr<Distribution1D> lightDistr;

    // Parameters light source sampling
    std::vector<int> nLightSamples;
    int nCausticPaths, nIndirectPaths;

    // Photon maps
    KdTree<Photon> *directMap;
    KdTree<Photon> *causticMap;
    KdTree<Photon> *indirectMap;
    KdTree<Photon> *volumeMap;
    KdTree<RadiancePhoton> *radianceMap;
};


VolPhotonIntegrator *CreateVolPhotonMapIntegrator(
        const ParamSet &params, std::shared_ptr<Sampler> sampler,
        std::shared_ptr<const Camera> camera);

}

#endif // PBRT_INTEGRATORS_VOLPHOTONMAP_H
