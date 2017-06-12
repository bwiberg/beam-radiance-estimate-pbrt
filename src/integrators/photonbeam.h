/**
 * Created by Benjamin Wiberg on 09/06/2017.
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PHOTONBEAM_H
#define PBRT_INTEGRATORS_PHOTONBEAM_H

// integrators/sppm.h*
#include "pbrt.h"
#include "integrator.h"
#include "camera.h"
#include "film.h"

namespace pbrt {

// SPPM Declarations
class PhotonBeamIntegrator : public Integrator {
public:
    // PhotonBeamIntegrator Public Methods
    PhotonBeamIntegrator(std::shared_ptr<const Camera> &camera, int nIterations,
                         int beamsPerIteraion, int maxDepth,
                         Float initialSearchRadius, Float alpha,
                         int writeFrequency,
                         bool rendersurfaces, bool rendermedia)
            : camera(camera),
              initialBeamRadius(initialSearchRadius),
              alpha(alpha),
              nIterations(nIterations),
              maxDepth(maxDepth),
              photonsPerIteration(beamsPerIteraion > 0
                                  ? beamsPerIteraion
                                  : camera->film->croppedPixelBounds.Area()),
              writeFrequency(writeFrequency),
              renderSurfaces(rendersurfaces),
              renderMedia(rendermedia)
    {}

    void Render(const Scene &scene);

private:
    // SPPMIntegrator Private Data
    std::shared_ptr<const Camera> camera;
    const Float initialBeamRadius;
    const Float alpha;
    const int nIterations;
    const int maxDepth;
    const int photonsPerIteration;
    const int writeFrequency;

    const bool renderSurfaces;
    const bool renderMedia;
};

Integrator *CreatePhotonBeamIntegrator(const ParamSet &params,
                                    std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_PHOTONBEAM_H
