
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


// photonbeambvh.cpp*
#include "photonbeambvh.h"
#include "paramset.h"
#include "stats.h"
#include "parallel.h"

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/BVH tree", treeBytes);
STAT_RATIO("PhotonBeamBVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);
STAT_COUNTER("PhotonBeamBVH/Interior nodes", interiorNodes);
STAT_COUNTER("PhotonBeamBVH/Leaf nodes", leafNodes);

// BVH Local Declarations
struct BVHBeamInfo {
    BVHBeamInfo() {}

    BVHBeamInfo(size_t primitiveNumber, const Bounds3f &bounds)
            : primitiveNumber(primitiveNumber),
              bounds(bounds),
              centroid(.5f * bounds.pMin + .5f * bounds.pMax) {}

    size_t primitiveNumber;
    Bounds3f bounds;
    Point3f centroid;
};

struct PBBVHBuildNode {
    // BVH2BuildNode Public Methods
    void InitLeaf(int first, int n, const Bounds3f &b) {
        firstPrimOffset = first;
        nPhotonBeams = n;
        bounds = b;
        children[0] = children[1] = nullptr;
        ++leafNodes;
        ++totalLeafNodes;
        totalPrimitives += n;
    }

    void InitInterior(int axis, PBBVHBuildNode *c0, PBBVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPhotonBeams = 0;
        ++interiorNodes;
    }

    Bounds3f bounds;
    PBBVHBuildNode *children[2];
    int splitAxis, firstPrimOffset, nPhotonBeams;
};

struct MortonPhotonBeam {
    int primitiveIndex;
    uint32_t mortonCode;
};

struct LBVH2Treelet {
    int startIndex, nPhotonBeams;
    PBBVHBuildNode *buildNodes;
};

struct LinearPBBVHNode {
    Bounds3f bounds;
    union {
        int photonBeamsOffset;   // leaf
        int secondChildOffset;  // interior
    };
    uint16_t nPhotonBeams;  // 0 -> interior node
    uint8_t axis;          // interior node: xyz
    uint8_t pad[1];        // ensure 32 byte total size
};

// BVH Utility Functions
inline uint32_t LeftShift3(uint32_t x) {
    CHECK_LE(x, (1 << 10));
    if (x == (1 << 10)) --x;
#ifdef PBRT_HAVE_BINARY_CONSTANTS
    x = (x | (x << 16)) & 0b00000011000000000000000011111111;
    // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x << 8)) & 0b00000011000000001111000000001111;
    // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x << 4)) & 0b00000011000011000011000011000011;
    // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x << 2)) & 0b00001001001001001001001001001001;
    // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#else
    x = (x | (x << 16)) & 0x30000ff;
    // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x << 8)) & 0x300f00f;
    // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x << 4)) & 0x30c30c3;
    // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x << 2)) & 0x9249249;
    // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#endif // PBRT_HAVE_BINARY_CONSTANTS
    return x;
}

inline uint32_t EncodeMorton3(const Vector3f &v) {
    CHECK_GE(v.x, 0);
    CHECK_GE(v.y, 0);
    CHECK_GE(v.z, 0);
    return (LeftShift3(v.z) << 2) | (LeftShift3(v.y) << 1) | LeftShift3(v.x);
}

static void RadixSort(std::vector<MortonPhotonBeam> *v) {
    std::vector<MortonPhotonBeam> tempVector(v->size());
    PBRT_CONSTEXPR int bitsPerPass = 6;
    PBRT_CONSTEXPR int nBits = 30;
    static_assert((nBits % bitsPerPass) == 0,
                  "Radix sort bitsPerPass must evenly divide nBits");
    PBRT_CONSTEXPR int nPasses = nBits / bitsPerPass;

    for (int pass = 0; pass < nPasses; ++pass) {
        // Perform one pass of radix sort, sorting _bitsPerPass_ bits
        int lowBit = pass * bitsPerPass;

        // Set in and out vector pointers for radix sort pass
        std::vector<MortonPhotonBeam> &in = (pass & 1) ? tempVector : *v;
        std::vector<MortonPhotonBeam> &out = (pass & 1) ? *v : tempVector;

        // Count number of zero bits in array for current radix sort bit
        PBRT_CONSTEXPR int nBuckets = 1 << bitsPerPass;
        int bucketCount[nBuckets] = {0};
        PBRT_CONSTEXPR int bitMask = (1 << bitsPerPass) - 1;
        for (const MortonPhotonBeam &mp : in) {
            int bucket = (mp.mortonCode >> lowBit) & bitMask;
            CHECK_GE(bucket, 0);
            CHECK_LT(bucket, nBuckets);
            ++bucketCount[bucket];
        }

        // Compute starting index in output array for each bucket
        int outIndex[nBuckets];
        outIndex[0] = 0;
        for (int i = 1; i < nBuckets; ++i)
            outIndex[i] = outIndex[i - 1] + bucketCount[i - 1];

        // Store sorted values in output array
        for (const MortonPhotonBeam &mp : in) {
            int bucket = (mp.mortonCode >> lowBit) & bitMask;
            out[outIndex[bucket]++] = mp;
        }
    }
    // Copy final result from _tempVector_, if needed
    if (nPasses & 1) std::swap(*v, tempVector);
}


//std::vector<std::shared_ptr<PhotonBeam>> Split(std::shared_ptr<PhotonBeam> &beam) {
//    Bounds3f bounds = beam->WorldBound()
//}


// BVH Method Definitions
PhotonBeamBVH::PhotonBeamBVH(const std::vector<std::shared_ptr<PhotonBeam>> &beams,
                             int maxPrimsInNode,
                             SplitMethod splitMethod)
        : maxPrimsInNode(std::min(255, maxPrimsInNode)),
          splitMethod(splitMethod),
          photonBeams(beams) {
    ProfilePhase _(Prof::PhotonBeamBVHConstruction);
    if (photonBeams.empty()) return;
    // Build BVH2 from _primitives_

    // Initialize _primitiveInfo_ array for primitives
    std::vector<BVHBeamInfo> primitiveInfo(photonBeams.size());
    for (size_t i = 0; i < photonBeams.size(); ++i)
        primitiveInfo[i] = {i, photonBeams[i]->WorldBound()};

    // Build BVH2 tree for primitives using _primitiveInfo_
    MemoryArena arena(1024 * 1024);
    int totalNodes = 0;
    std::vector<std::shared_ptr<PhotonBeam>> orderedBeams;
    orderedBeams.reserve(photonBeams.size());
    PBBVHBuildNode *root;
    if (splitMethod == SplitMethod::HLBVH)
        root = HLBVH2Build(arena, primitiveInfo, &totalNodes, orderedBeams);
    else
        root = recursiveBuild(arena, primitiveInfo, 0, photonBeams.size(),
                              &totalNodes, orderedBeams);
    photonBeams.swap(orderedBeams);
    LOG(INFO) << StringPrintf("BVH2 created with %d nodes for %d "
                                      "primitives (%.2f MB)", totalNodes,
                              (int) photonBeams.size(),
                              float(totalNodes * sizeof(LinearPBBVHNode)) /
                              (1024.f * 1024.f));

    // Compute representation of depth-first traversal of BVH2 tree
    treeBytes += totalNodes * sizeof(LinearPBBVHNode) + sizeof(*this) +
                 photonBeams.size() * sizeof(photonBeams[0]);
    nodes = AllocAligned<LinearPBBVHNode>(totalNodes);
    int offset = 0;
    flattenBVH2Tree(root, &offset);
    CHECK_EQ(totalNodes, offset);
}

Bounds3f PhotonBeamBVH::WorldBound() const {
    return nodes ? nodes[0].bounds : Bounds3f();
}

struct BucketInfo {
    int count = 0;
    Bounds3f bounds;
};

PBBVHBuildNode *PhotonBeamBVH::recursiveBuild(
        MemoryArena &arena, std::vector<BVHBeamInfo> &primitiveInfo, int start,
        int end, int *totalNodes,
        std::vector<std::shared_ptr<PhotonBeam>> &orderedBeams) {
    CHECK_NE(start, end);
    PBBVHBuildNode *node = arena.Alloc<PBBVHBuildNode>();
    (*totalNodes)++;
    // Compute bounds of all primitives in BVH2 node
    Bounds3f bounds;
    for (int i = start; i < end; ++i)
        bounds = Union(bounds, primitiveInfo[i].bounds);
    int nPhotonBeams = end - start;
    if (nPhotonBeams == 1) {
        // Create leaf _BVH2BuildNode_
        int firstPrimOffset = orderedBeams.size();
        for (int i = start; i < end; ++i) {
            int primNum = primitiveInfo[i].primitiveNumber;
            orderedBeams.push_back(photonBeams[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPhotonBeams, bounds);
        return node;
    } else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        Bounds3f centroidBounds;
        for (int i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        int mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // Create leaf _BVH2BuildNode_
            int firstPrimOffset = orderedBeams.size();
            for (int i = start; i < end; ++i) {
                int primNum = primitiveInfo[i].primitiveNumber;
                orderedBeams.push_back(photonBeams[primNum]);
            }
            node->InitLeaf(firstPrimOffset, nPhotonBeams, bounds);
            return node;
        } else {
            // Partition primitives based on _splitMethod_
            switch (splitMethod) {
                case SplitMethod::Middle: {
                    // Partition primitives through node's midpoint
                    Float pmid =
                            (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
                    BVHBeamInfo *midPtr = std::partition(
                            &primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                            [dim, pmid](const BVHBeamInfo &pi) {
                                return pi.centroid[dim] < pmid;
                            });
                    mid = midPtr - &primitiveInfo[0];
                    // For lots of prims with large overlapping bounding boxes, this
                    // may fail to partition; in that case don't break and fall
                    // through
                    // to EqualCounts.
                    if (mid != start && mid != end) break;
                }
                case SplitMethod::EqualCounts: {
                    // Partition primitives into equally-sized subsets
                    mid = (start + end) / 2;
                    std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
                                     &primitiveInfo[end - 1] + 1,
                                     [dim](const BVHBeamInfo &a,
                                           const BVHBeamInfo &b) {
                                         return a.centroid[dim] < b.centroid[dim];
                                     });
                    break;
                }
                case SplitMethod::SAH:
                default: {
                    // Partition primitives using approximate SAH
                    if (nPhotonBeams <= 2) {
                        // Partition primitives into equally-sized subsets
                        mid = (start + end) / 2;
                        std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
                                         &primitiveInfo[end - 1] + 1,
                                         [dim](const BVHBeamInfo &a,
                                               const BVHBeamInfo &b) {
                                             return a.centroid[dim] <
                                                    b.centroid[dim];
                                         });
                    } else {
                        // Allocate _BucketInfo_ for SAH partition buckets
                        PBRT_CONSTEXPR int nBuckets = 12;
                        BucketInfo buckets[nBuckets];

                        // Initialize _BucketInfo_ for SAH partition buckets
                        for (int i = start; i < end; ++i) {
                            int b = nBuckets *
                                    centroidBounds.Offset(
                                            primitiveInfo[i].centroid)[dim];
                            if (b == nBuckets) b = nBuckets - 1;
                            CHECK_GE(b, 0);
                            CHECK_LT(b, nBuckets);
                            buckets[b].count++;
                            buckets[b].bounds =
                                    Union(buckets[b].bounds, primitiveInfo[i].bounds);
                        }

                        // Compute costs for splitting after each bucket
                        Float cost[nBuckets - 1];
                        for (int i = 0; i < nBuckets - 1; ++i) {
                            Bounds3f b0, b1;
                            int count0 = 0, count1 = 0;
                            for (int j = 0; j <= i; ++j) {
                                b0 = Union(b0, buckets[j].bounds);
                                count0 += buckets[j].count;
                            }
                            for (int j = i + 1; j < nBuckets; ++j) {
                                b1 = Union(b1, buckets[j].bounds);
                                count1 += buckets[j].count;
                            }
                            cost[i] = 1 +
                                      (count0 * b0.SurfaceArea() +
                                       count1 * b1.SurfaceArea()) /
                                      bounds.SurfaceArea();
                        }

                        // Find bucket to split at that minimizes SAH metric
                        Float minCost = cost[0];
                        int minCostSplitBucket = 0;
                        for (int i = 1; i < nBuckets - 1; ++i) {
                            if (cost[i] < minCost) {
                                minCost = cost[i];
                                minCostSplitBucket = i;
                            }
                        }

                        // Either create leaf or split primitives at selected SAH
                        // bucket
                        Float leafCost = nPhotonBeams;
                        if (nPhotonBeams > maxPrimsInNode || minCost < leafCost) {
                            BVHBeamInfo *pmid = std::partition(
                                    &primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                                    [=](const BVHBeamInfo &pi) {
                                        int b = nBuckets *
                                                centroidBounds.Offset(pi.centroid)[dim];
                                        if (b == nBuckets) b = nBuckets - 1;
                                        CHECK_GE(b, 0);
                                        CHECK_LT(b, nBuckets);
                                        return b <= minCostSplitBucket;
                                    });
                            mid = pmid - &primitiveInfo[0];
                        } else {
                            // Create leaf _BVH2BuildNode_
                            int firstPrimOffset = orderedBeams.size();
                            for (int i = start; i < end; ++i) {
                                int primNum = primitiveInfo[i].primitiveNumber;
                                orderedBeams.push_back(photonBeams[primNum]);
                            }
                            node->InitLeaf(firstPrimOffset, nPhotonBeams, bounds);
                            return node;
                        }
                    }
                    break;
                }
            }
            node->InitInterior(dim,
                               recursiveBuild(arena, primitiveInfo, start, mid,
                                              totalNodes, orderedBeams),
                               recursiveBuild(arena, primitiveInfo, mid, end,
                                              totalNodes, orderedBeams));
        }
    }
    return node;
}

PBBVHBuildNode *PhotonBeamBVH::HLBVH2Build(
        MemoryArena &arena, const std::vector<BVHBeamInfo> &primitiveInfo,
        int *totalNodes,
        std::vector<std::shared_ptr<PhotonBeam>> &orderedBeams) const {
    // Compute bounding box of all primitive centroids
    Bounds3f bounds;
    for (const BVHBeamInfo &pi : primitiveInfo)
        bounds = Union(bounds, pi.centroid);

    // Compute Morton indices of primitives
    std::vector<MortonPhotonBeam> mortonPrims(primitiveInfo.size());
    ParallelFor([&](int i) {
        // Initialize _mortonPrims[i]_ for _i_th primitive
        PBRT_CONSTEXPR int mortonBits = 10;
        PBRT_CONSTEXPR int mortonScale = 1 << mortonBits;
        mortonPrims[i].primitiveIndex = primitiveInfo[i].primitiveNumber;
        Vector3f centroidOffset = bounds.Offset(primitiveInfo[i].centroid);
        mortonPrims[i].mortonCode = EncodeMorton3(centroidOffset * mortonScale);
    }, primitiveInfo.size(), 512);

    // Radix sort primitive Morton indices
    RadixSort(&mortonPrims);

    // Create LBVH2 treelets at bottom of BVH2

    // Find intervals of primitives for each treelet
    std::vector<LBVH2Treelet> treeletsToBuild;
    for (int start = 0, end = 1; end <= (int) mortonPrims.size(); ++end) {
#ifdef PBRT_HAVE_BINARY_CONSTANTS
        uint32_t mask = 0b00111111111111000000000000000000;
#else
        uint32_t mask = 0x3ffc0000;
#endif
        if (end == (int) mortonPrims.size() ||
            ((mortonPrims[start].mortonCode & mask) !=
             (mortonPrims[end].mortonCode & mask))) {
            // Add entry to _treeletsToBuild_ for this treelet
            int nPhotonBeams = end - start;
            int maxBVH2Nodes = 2 * nPhotonBeams;
            PBBVHBuildNode *nodes = arena.Alloc<PBBVHBuildNode>(maxBVH2Nodes, false);
            treeletsToBuild.push_back({start, nPhotonBeams, nodes});
            start = end;
        }
    }

    // Create LBVH2s for treelets in parallel
    std::atomic<int> atomicTotal(0), orderedBeamsOffset(0);
    orderedBeams.resize(photonBeams.size());
    ParallelFor([&](int i) {
        // Generate _i_th LBVH2 treelet
        int nodesCreated = 0;
        const int firstBitIndex = 29 - 12;
        LBVH2Treelet &tr = treeletsToBuild[i];
        tr.buildNodes =
                emitLBVH2(tr.buildNodes, primitiveInfo, &mortonPrims[tr.startIndex],
                          tr.nPhotonBeams, &nodesCreated, orderedBeams,
                          &orderedBeamsOffset, firstBitIndex);
        atomicTotal += nodesCreated;
    }, treeletsToBuild.size());
    *totalNodes = atomicTotal;

    // Create and return SAH BVH2 from LBVH2 treelets
    std::vector<PBBVHBuildNode *> finishedTreelets;
    finishedTreelets.reserve(treeletsToBuild.size());
    for (LBVH2Treelet &treelet : treeletsToBuild)
        finishedTreelets.push_back(treelet.buildNodes);
    return buildUpperSAH(arena, finishedTreelets, 0, finishedTreelets.size(),
                         totalNodes);
}

PBBVHBuildNode *PhotonBeamBVH::emitLBVH2(
        PBBVHBuildNode *&buildNodes,
        const std::vector<BVHBeamInfo> &primitiveInfo,
        MortonPhotonBeam *mortonPrims, int nPhotonBeams, int *totalNodes,
        std::vector<std::shared_ptr<PhotonBeam>> &orderedBeams,
        std::atomic<int> *orderedBeamsOffset, int bitIndex) const {
    CHECK_GT(nPhotonBeams, 0);
    if (bitIndex == -1 || nPhotonBeams < maxPrimsInNode) {
        // Create and return leaf node of LBVH2 treelet
        (*totalNodes)++;
        PBBVHBuildNode *node = buildNodes++;
        Bounds3f bounds;
        int firstPrimOffset = orderedBeamsOffset->fetch_add(nPhotonBeams);
        for (int i = 0; i < nPhotonBeams; ++i) {
            int primitiveIndex = mortonPrims[i].primitiveIndex;
            orderedBeams[firstPrimOffset + i] = photonBeams[primitiveIndex];
            bounds = Union(bounds, primitiveInfo[primitiveIndex].bounds);
        }
        node->InitLeaf(firstPrimOffset, nPhotonBeams, bounds);
        return node;
    } else {
        int mask = 1 << bitIndex;
        // Advance to next subtree level if there's no LBVH2 split for this bit
        if ((mortonPrims[0].mortonCode & mask) ==
            (mortonPrims[nPhotonBeams - 1].mortonCode & mask))
            return emitLBVH2(buildNodes, primitiveInfo, mortonPrims, nPhotonBeams,
                             totalNodes, orderedBeams, orderedBeamsOffset,
                             bitIndex - 1);

        // Find LBVH2 split point for this dimension
        int searchStart = 0, searchEnd = nPhotonBeams - 1;
        while (searchStart + 1 != searchEnd) {
            CHECK_NE(searchStart, searchEnd);
            int mid = (searchStart + searchEnd) / 2;
            if ((mortonPrims[searchStart].mortonCode & mask) ==
                (mortonPrims[mid].mortonCode & mask))
                searchStart = mid;
            else {
                CHECK_EQ(mortonPrims[mid].mortonCode & mask,
                         mortonPrims[searchEnd].mortonCode & mask);
                searchEnd = mid;
            }
        }
        int splitOffset = searchEnd;
        CHECK_LE(splitOffset, nPhotonBeams - 1);
        CHECK_NE(mortonPrims[splitOffset - 1].mortonCode & mask,
                 mortonPrims[splitOffset].mortonCode & mask);

        // Create and return interior LBVH2 node
        (*totalNodes)++;
        PBBVHBuildNode *node = buildNodes++;
        PBBVHBuildNode *lbvh[2] = {
                emitLBVH2(buildNodes, primitiveInfo, mortonPrims, splitOffset,
                          totalNodes, orderedBeams, orderedBeamsOffset,
                          bitIndex - 1),
                emitLBVH2(buildNodes, primitiveInfo, &mortonPrims[splitOffset],
                          nPhotonBeams - splitOffset, totalNodes, orderedBeams,
                          orderedBeamsOffset, bitIndex - 1)};
        int axis = bitIndex % 3;
        node->InitInterior(axis, lbvh[0], lbvh[1]);
        return node;
    }
}

PBBVHBuildNode *PhotonBeamBVH::buildUpperSAH(MemoryArena &arena,
                                             std::vector<PBBVHBuildNode *> &treeletRoots,
                                             int start, int end,
                                             int *totalNodes) const {
    CHECK_LT(start, end);
    int nNodes = end - start;
    if (nNodes == 1) return treeletRoots[start];
    (*totalNodes)++;
    PBBVHBuildNode *node = arena.Alloc<PBBVHBuildNode>();

    // Compute bounds of all nodes under this HLBVH2 node
    Bounds3f bounds;
    for (int i = start; i < end; ++i)
        bounds = Union(bounds, treeletRoots[i]->bounds);

    // Compute bound of HLBVH2 node centroids, choose split dimension _dim_
    Bounds3f centroidBounds;
    for (int i = start; i < end; ++i) {
        Point3f centroid =
                (treeletRoots[i]->bounds.pMin + treeletRoots[i]->bounds.pMax) *
                0.5f;
        centroidBounds = Union(centroidBounds, centroid);
    }
    int dim = centroidBounds.MaximumExtent();
    // FIXME: if this hits, what do we need to do?
    // Make sure the SAH split below does something... ?
    CHECK_NE(centroidBounds.pMax[dim], centroidBounds.pMin[dim]);

    // Allocate _BucketInfo_ for SAH partition buckets
    PBRT_CONSTEXPR int nBuckets = 12;
    struct BucketInfo {
        int count = 0;
        Bounds3f bounds;
    };
    BucketInfo buckets[nBuckets];

    // Initialize _BucketInfo_ for HLBVH2 SAH partition buckets
    for (int i = start; i < end; ++i) {
        Float centroid = (treeletRoots[i]->bounds.pMin[dim] +
                          treeletRoots[i]->bounds.pMax[dim]) *
                         0.5f;
        int b =
                nBuckets * ((centroid - centroidBounds.pMin[dim]) /
                            (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
        if (b == nBuckets) b = nBuckets - 1;
        CHECK_GE(b, 0);
        CHECK_LT(b, nBuckets);
        buckets[b].count++;
        buckets[b].bounds = Union(buckets[b].bounds, treeletRoots[i]->bounds);
    }

    // Compute costs for splitting after each bucket
    Float cost[nBuckets - 1];
    for (int i = 0; i < nBuckets - 1; ++i) {
        Bounds3f b0, b1;
        int count0 = 0, count1 = 0;
        for (int j = 0; j <= i; ++j) {
            b0 = Union(b0, buckets[j].bounds);
            count0 += buckets[j].count;
        }
        for (int j = i + 1; j < nBuckets; ++j) {
            b1 = Union(b1, buckets[j].bounds);
            count1 += buckets[j].count;
        }
        cost[i] = .125f +
                  (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) /
                  bounds.SurfaceArea();
    }

    // Find bucket to split at that minimizes SAH metric
    Float minCost = cost[0];
    int minCostSplitBucket = 0;
    for (int i = 1; i < nBuckets - 1; ++i) {
        if (cost[i] < minCost) {
            minCost = cost[i];
            minCostSplitBucket = i;
        }
    }

    // Split nodes and create interior HLBVH2 SAH node
    PBBVHBuildNode **pmid = std::partition(
            &treeletRoots[start], &treeletRoots[end - 1] + 1,
            [=](const PBBVHBuildNode *node) {
                Float centroid =
                        (node->bounds.pMin[dim] + node->bounds.pMax[dim]) * 0.5f;
                int b = nBuckets *
                        ((centroid - centroidBounds.pMin[dim]) /
                         (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
                if (b == nBuckets) b = nBuckets - 1;
                CHECK_GE(b, 0);
                CHECK_LT(b, nBuckets);
                return b <= minCostSplitBucket;
            });
    int mid = pmid - &treeletRoots[0];
    CHECK_GT(mid, start);
    CHECK_LT(mid, end);
    node->InitInterior(
            dim, this->buildUpperSAH(arena, treeletRoots, start, mid, totalNodes),
            this->buildUpperSAH(arena, treeletRoots, mid, end, totalNodes));
    return node;
}

int PhotonBeamBVH::flattenBVH2Tree(PBBVHBuildNode *node, int *offset) {
    LinearPBBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    int myOffset = (*offset)++;
    if (node->nPhotonBeams > 0) {
        CHECK(!node->children[0] && !node->children[1]);
        CHECK_LT(node->nPhotonBeams, 65536);
        linearNode->photonBeamsOffset = node->firstPrimOffset;
        linearNode->nPhotonBeams = node->nPhotonBeams;
    } else {
        // Create interior flattened BVH2 node
        linearNode->axis = node->splitAxis;
        linearNode->nPhotonBeams = 0;
        flattenBVH2Tree(node->children[0], offset);
        linearNode->secondChildOffset =
                flattenBVH2Tree(node->children[1], offset);
    }
    return myOffset;
}

PhotonBeamBVH::~PhotonBeamBVH() { FreeAligned(nodes); }

std::vector<std::shared_ptr<PhotonBeam>> PhotonBeamBVH::Intersect(const Ray &ray) const {
    std::vector<std::shared_ptr<PhotonBeam>> beams;
    if (!nodes) return beams;

    ProfilePhase p(Prof::PhotonBeamBVHQuery);
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};

    // Follow ray through BVH2 nodes to find primitive intersections
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[64];
    while (true) {
        const LinearPBBVHNode *node = &nodes[currentNodeIndex];
        // Check ray against BVH2 node
        if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
            if (node->nPhotonBeams > 0) {
                // Intersect ray with primitives in leaf BVH2 node
                for (int i = 0; i < node->nPhotonBeams; ++i)
                    beams.push_back(photonBeams[node->photonBeamsOffset + i]);
                if (toVisitOffset == 0) break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                // Put far BVH2 node on _nodesToVisit_ stack, advance to near
                // node
                if (dirIsNeg[node->axis]) {
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
        } else {
            if (toVisitOffset == 0) break;
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }
    }
    return beams;
}
}  // namespace pbrt
