
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_BVH2_H
#define PBRT_ACCELERATORS_BVH2_H

// accelerators/bvh.h*
#include "pbrt.h"
#include "geometry.h"
#include "spectrum.h"
#include <atomic>

namespace pbrt {
struct PhotonBeam {
    PhotonBeam(Point3f const &start, Point3f const &end,
               Float radius, Spectrum const &powerstart, Spectrum const &powerend)
            : start(start), end(end),
              radius(radius),
              powerStart(powerstart),
              powerEnd(powerend){}

    Point3f start, end;
    Float radius;
    Spectrum powerStart, powerEnd;

    Bounds3f WorldBound() const {
        Vector3f dir = end - start;
        const Point3f center = start + dir / 2;
        const Float len = dir.Length();
        dir /= len;

        Vector3f size;
        size.x = dir.x * len + 2 * radius * sqrt(1 - dir.x * dir.x);
        size.y = dir.y * len + 2 * radius * sqrt(1 - dir.y * dir.y);
        size.z = dir.z * len + 2 * radius * sqrt(1 - dir.z * dir.z);

        return Bounds3f(center - size / 2, center + size / 2);
    }
};

struct PBBVHBuildNode;

// BVH Forward Declarations
struct BVHBeamInfo;
struct MortonPhotonBeam;
struct LinearPBBVHNode;

// BVH Declarations
class PhotonBeamBVH {
public:
    // BVH Public Types
    enum class SplitMethod {
        SAH, HLBVH, Middle, EqualCounts
    };

    // BVH Public Methods
    PhotonBeamBVH(std::vector<std::shared_ptr<PhotonBeam>> && beams,
                  int maxPrimsInNode = 1,
                  SplitMethod splitMethod = SplitMethod::SAH);

    Bounds3f WorldBound() const;

    ~PhotonBeamBVH();

    std::vector<std::shared_ptr<PhotonBeam>> Intersect(const Ray &ray) const;

private:
    // BVH Private Methods
    PBBVHBuildNode *recursiveBuild(
            MemoryArena &arena, std::vector<BVHBeamInfo> &photonBeams,
            int start, int end, int *totalNodes,
            std::vector<std::shared_ptr<PhotonBeam>> &orderedPrims);

    PBBVHBuildNode *HLBVH2Build(
            MemoryArena &arena, const std::vector<BVHBeamInfo> &photonBeams,
            int *totalNodes,
            std::vector<std::shared_ptr<PhotonBeam>> &orderedPrims) const;

    PBBVHBuildNode *emitLBVH2(
            PBBVHBuildNode *&buildNodes,
            const std::vector<BVHBeamInfo> &primitiveInfo,
            MortonPhotonBeam *mortonPrims, int nPrimitives, int *totalNodes,
            std::vector<std::shared_ptr<PhotonBeam>> &orderedPrims,
            std::atomic<int> *orderedPrimsOffset, int bitIndex) const;

    PBBVHBuildNode *buildUpperSAH(MemoryArena &arena,
                                  std::vector<PBBVHBuildNode *> &treeletRoots,
                                  int start, int end, int *totalNodes) const;

    int flattenBVH2Tree(PBBVHBuildNode *node, int *offset);

    // BVH Private Data
    const int maxPrimsInNode;
    const SplitMethod splitMethod;
    std::vector<std::shared_ptr<PhotonBeam>> photonBeams;
    LinearPBBVHNode *nodes = nullptr;
};

std::shared_ptr<PhotonBeamBVH> CreateBVHerator(
        const std::vector<std::shared_ptr<PhotonBeam>> &prims, const ParamSet &ps);

}  // namespace pbrt

#endif  // PBRT_ACCELERATORS_BVH2_H
