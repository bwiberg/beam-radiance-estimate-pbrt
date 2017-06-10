/**
 * Created by Benjamin Wiberg on 09/06/2017.
 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_VSPPM_H
#define PBRT_INTEGRATORS_VSPPM_H

// integrators/sppm.h*
#include "pbrt.h"
#include "integrator.h"
#include "camera.h"
#include "film.h"

namespace pbrt {

enum class PhotonType {
    POINT = 0,
    BEAM = 1
};

// SPPM Declarations
class VolSPPMIntegrator : public Integrator {
public:
    // VolSPPMIntegrator Public Methods
    VolSPPMIntegrator(std::shared_ptr<const Camera> &camera, int nIterations,
                      int photonsPerIteration, int maxDepth,
                      Float initialSearchRadius, int writeFrequency,
                      PhotonType photonType)
            : camera(camera),
              initialSearchRadius(initialSearchRadius),
              nIterations(nIterations),
              maxDepth(maxDepth),
              photonsPerIteration(photonsPerIteration > 0
                                  ? photonsPerIteration
                                  : camera->film->croppedPixelBounds.Area()),
              writeFrequency(writeFrequency),
              photonType(photonType){}

    void Render(const Scene &scene);

private:
    // SPPMIntegrator Private Data
    std::shared_ptr<const Camera> camera;
    const Float initialSearchRadius;
    const int nIterations;
    const int maxDepth;
    const int photonsPerIteration;
    const int writeFrequency;

    const PhotonType photonType;
};

Integrator *CreateVolSPPMIntegrator(const ParamSet &params,
                                    std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_VSPPM_H
