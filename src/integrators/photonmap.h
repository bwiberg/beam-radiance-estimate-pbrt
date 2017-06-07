
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

#ifndef PBRT_INTEGRATORS_PHOTONMAP_H
#define PBRT_INTEGRATORS_PHOTONMAP_H

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

// PhotonIntegrator Declarations
class PhotonIntegrator : public SamplerIntegrator {
public:
    // PhotonIntegrator Public Methods
    PhotonIntegrator(std::shared_ptr<const Camera> camera,
                     std::shared_ptr<Sampler> sampler,
                     const Bounds2i &pixelBounds,
                     int ncaus, int nindir, bool reqphotons,
                     int nLookup, int maxspecdepth,
                     int maxphotondepth, float maxdist, bool finalGather, int gatherSamples,
                     float ga);

    ~PhotonIntegrator();

    Spectrum Li(const RayDifferential &ray,
                const Scene &scene,
                Sampler &sampler,
                MemoryArena &arena,
                int depth = 0) const;

    void Preprocess(const Scene &scene, Sampler &sampler);

private:
    void ComputeLightSamples(Scene const& scene, Sampler &sampler);

    void ComputeLightPhotonDistrib(Scene const& scene);

    void ShootPhotons(Scene const& scene,
                      std::vector<Photon> &causticPhotons,
                      std::vector<Photon> &directPhotons,
                      std::vector<Photon> &indirectPhotons,
                      std::vector<RadiancePhoton> &radiancePhotons,
                      std::vector<Spectrum> &rpReflectances, std::vector<Spectrum> &rpTransmittances,
                      int &nDirectPaths, uint32_t &nShot,
                      ProgressReporter &progress);

    KdTree<Photon> *BuildPhotonMaps(std::vector<Photon> const& directPhotons,
                                    std::vector<Photon> const& causticPhotons,
                                    std::vector<Photon> const& indirectPhotons);

    // PhotonIntegrator Private Methods
    void ComputePhotonRadiances(std::vector<RadiancePhoton> &radiancePhotons,
                                std::vector<Spectrum> const& rpReflectances,
                                std::vector<Spectrum> const& rpTransmittances,
                                KdTree<Photon> const* directMap,
                                int nDirectPaths,
                                ProgressReporter &progress);

    // PhotonIntegrator Private Data
    uint32_t nCausticPhotonsWanted, nIndirectPhotonsWanted, nLookup;
    bool requirePhotons;
    Float maxDistSquared;
    int maxSpecularDepth, maxPhotonDepth;
    bool finalGather;
    int gatherSamples;
    Float cosGatherAngle;

    const std::string lightSampleStrategy;
    std::unique_ptr<Distribution1D> lightDistr;

    // Declare sample parameters for light source sampling
    std::vector<int> nLightSamples;
    int nCausticPaths, nIndirectPaths;
    KdTree<Photon> *causticMap;
    KdTree<Photon> *indirectMap;
    KdTree<RadiancePhoton> *radianceMap;
};


PhotonIntegrator *CreatePhotonMapIntegrator(
        const ParamSet &params, std::shared_ptr<Sampler> sampler,
        std::shared_ptr<const Camera> camera);

}

#endif // PBRT_INTEGRATORS_PHOTONMAP_H
