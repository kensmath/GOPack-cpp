#include "gopack/Packer.h"

#include <algorithm>
#include <cstdio>
#include <deque>
#include <random>
#include <stdexcept>

namespace gopack {

void Packer::cleanse() {
    fileName = "noname";
    alpha = 0;
    gamma = 0;
    hes = Geometry::Euclidean;
    nodeCount = 0;
    intCount = 0;
    bdryCount = 0;
    orphanCount = 0;
    bdryFlags.clear();
    flowers.clear();
    vNum.clear();
    vAims.clear();
    vlist.clear();
    edgeCount = 0;
    faceCount = 0;

    intVerts.clear();
    bdryList.clear();
    orphanVerts.clear();
    layoutVerts.clear();

    origRadii.clear();
    origCenters.clear();
    localradii.clear();
    localcenters.clear();
    radii.clear();
    centers.clear();
    v2indx.clear();
    indx2v.clear();

    corners.clear();
    sides.clear();
    mode = 1;
    angsumMonitor.clear();
    l2Monitor.clear();
    visErrMonitor.clear();
    ticMonitor.clear();
}

Index Packer::farVert(const std::vector<Index>& seeds) const {
    if (seeds.empty()) {
        return -1;
    }

    std::vector<int> marks(static_cast<size_t>(nodeCount) + 1, 0);
    std::vector<Index> nextlist = seeds;
    int gennum = 1;
    Index farvertResult = seeds[0];
    Index nothits = nodeCount;

    while (!nextlist.empty() && (gennum < 10 || nothits < 100)) {
        std::vector<Index> curr;
        for (Index v : nextlist) {
            if (v > 0) curr.push_back(v);
        }
        std::vector<Index> next;
        for (Index k : curr) {
            if (marks[k] < 1) {
                nothits--;
            }
            marks[k] = gennum;
            farvertResult = k;
            for (Index p : flowers[k]) {
                if (marks[p] == 0) {
                    next.push_back(p);
                }
            }
        }
        nextlist = std::move(next);
        gennum++;
    }

    if (nothits <= 100) {
        return farvertResult;
    }

    // Choose randomly among unreached vertices.
    static thread_local std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<Index> dist(1, nodeCount);
    while (true) {
        Index k = dist(rng);
        if (marks[k] == 0) {
            return k;
        }
    }
}

Index Packer::complexCount() {
    // determine alpha/gamma validity
    Index a = alpha;
    if (a < 1 || a > nodeCount || flowers[a].front() != flowers[a].back()) {
        a = -1;
    }
    Index g = gamma;
    if (g == a) {
        g = -1;
    }

    bool hasbdry = false;
    std::vector<int> tmpBdryFlags(static_cast<size_t>(nodeCount) + 1, 0);
    for (Index v = 1; v <= nodeCount; ++v) {
        if (flowers[v].front() == flowers[v].back()) {
            tmpBdryFlags[v] = 0;
            if (a == -1) a = v;
        } else {
            tmpBdryFlags[v] = 1;
            if (g == -1) g = v;
            hasbdry = true;
        }
    }

    if (a == -1) {
        std::fprintf(stderr, "Error: complex has no interior vertex\n");
        return -1;
    }
    alpha = a;

    bdryList.clear();
    layoutVerts.clear();
    orphanVerts.clear();
    orphanEdges.clear();

    if (!hasbdry) {
        // ---- spherical case ----
        hes = Geometry::Spherical;
        Index av = farVert({alpha});
        const auto& flower = flowers[av];
        Index b = flower[1]; // flower(2) in 1-indexed MATLAB -> index 1 here (0-indexed within the vector, but vertex labels are still 1-indexed values)
        Index c = flower[0];
        bdryList = {av, b, c, av};
        gamma = av;
        intVerts.clear();
        for (Index v = 1; v <= nodeCount; ++v) {
            if (v != av && v != b && v != c) {
                intVerts.push_back(v);
            }
        }
        orphanVerts.clear();
    } else {
        // ---- non-spherical: BFS from alpha over interior-only edges ----
        std::vector<int> status(static_cast<size_t>(nodeCount) + 1, 0);
        std::deque<Index> hitlist;
        hitlist.push_back(alpha);
        status[alpha] = -1;
        intVerts.clear();
        intVerts.push_back(alpha);

        while (!hitlist.empty()) {
            Index v = hitlist.front();
            hitlist.pop_front();
            for (Index w : flowers[v]) {
                if (tmpBdryFlags[w] == 0 && status[w] == 0) {
                    status[w] = -1;
                    hitlist.push_back(w);
                    intVerts.push_back(w);
                }
            }
            status[v] = 1;
        }

        std::vector<Index> tmpBdryV;
        for (Index intV : intVerts) {
            for (Index w : flowers[intV]) {
                if (status[w] == 0) {
                    status[w] = -1;
                    tmpBdryV.push_back(w);
                }
            }
        }

        if (tmpBdryV.empty()) {
            std::fprintf(stderr, "Error forming bdryList: no boundary vertices found\n");
            return -1;
        }

        Index firstbdry = gamma;
        if (firstbdry < 1 || firstbdry > nodeCount || status[firstbdry] != -1) {
            firstbdry = tmpBdryV[0];
            gamma = firstbdry;
        }

        bdryList.clear();
        bdryList.push_back(firstbdry);

        Index nextb = 0;
        for (Index w : flowers[firstbdry]) {
            if (status[w] == -1) {
                bdryList.push_back(w);
                nextb = w;
                break;
            }
        }

        const size_t maxTick = 2 * tmpBdryV.size();
        while (nextb != firstbdry && bdryList.size() < maxTick) {
            Index found = 0;
            for (Index w : flowers[nextb]) {
                if (status[w] == -1) {
                    bdryList.push_back(w);
                    found = w;
                    break;
                }
            }
            nextb = found;
            if (found == 0) break; // no progress; avoid infinite loop
        }
        // Bug fix vs. complex_count.m: the source prints these two
        // diagnostics but has no early return here, so it falls straight
        // through and returns a "successful" nodeCount even though bdryList
        // is now known to be an incomplete/malformed loop (didn't close
        // back on firstbdry, or ran away past every genuine boundary
        // vertex). Every caller further down the line -- pruneComplex(),
        // indxMatrices(), and eventually layoutCenters()'s Eigen sparse
        // assembly -- assumes bdryList is a complete, correctly-closed
        // cycle; feeding them a partial one produces vertices that are
        // real boundary neighbors of the interior (they're in tmpBdryV) but
        // never made it into bdryList/rimVerts, so indxMatrices() later
        // treats them as having v2indx==0 -- silently miscategorized as an
        // "interior" neighbor with column index 0, corrupting the sparse
        // triplet lists layoutCenters() feeds straight to Eigen. In a
        // Release build that reads as an unpredictable crash deep inside
        // riffle(), not a clean, diagnosable error. This function's own doc
        // comment in Packer.h already documents a "-1 on error" contract
        // (matching the "no interior vertex" case above) -- honor it here
        // too instead of silently returning a corrupt complex.
        if (nextb != firstbdry || bdryList.back() != bdryList.front()) {
            std::fprintf(stderr, "Error forming bdryList\n");
            return -1;
        }
        if (bdryList.size() >= maxTick) {
            std::fprintf(stderr, "Error in bdryList, too long\n");
            return -1;
        }

        orphanVerts.clear();
        for (Index i = 1; i <= nodeCount; ++i) {
            if (status[i] == 0) {
                orphanVerts.push_back(i);
            }
        }
    }

    intCount = static_cast<Index>(intVerts.size());
    bdryCount = static_cast<Index>(bdryList.size()) - 1;
    orphanCount = static_cast<Index>(orphanVerts.size());

    edgeCount = 0;
    faceCount = 0;
    vNum.assign(static_cast<size_t>(nodeCount) + 1, 0);
    bdryFlags.assign(static_cast<size_t>(nodeCount) + 1, 0);
    Index totNum = 0;
    for (Index v = 1; v <= nodeCount; ++v) {
        const auto& flower = flowers[v];
        Index num = static_cast<Index>(flower.size()) - 1;
        vNum[v] = num;
        totNum += num;
        for (Index k = 0; k < num; ++k) {
            Index w = flower[k];
            if (w > v) edgeCount++;
        }
        if (flower.front() != flower.back()) {
            Index w = flower.back();
            if (w > v) edgeCount++;
            bdryFlags[w] = 1;
        }
    }
    faceCount = totNum / 3;

    // ---- orphanEdges ----
    orphanEdges.clear();
    if (!orphanVerts.empty()) {
        orphanEdges.assign(static_cast<size_t>(orphanCount), {0, 0});
        std::vector<Index> ctlg(static_cast<size_t>(nodeCount) + 1, 0);
        for (Index v : intVerts) ctlg[v] = -2;
        for (size_t j = 0; j + 1 < bdryList.size(); ++j) ctlg[bdryList[j]] = -1;

        int tick = 0;
        for (size_t j = 0; j < orphanVerts.size(); ++j) {
            Index v = orphanVerts[j];
            const auto& flower = flowers[v];
            for (Index k = 0; k + 1 < static_cast<Index>(flower.size()); ++k) {
                Index m = flower[k];
                Index n = flower[k + 1];
                if (ctlg[m] == -1 && ctlg[n] == -1) {
                    orphanEdges[j] = {n, m};
                    ctlg[v] = static_cast<Index>(j) + 1;
                    tick++;
                    break;
                }
            }
        }

        int hit = 1;
        while (hit > 0) {
            hit = 0;
            for (size_t j = 0; j < orphanVerts.size(); ++j) {
                Index v = orphanVerts[j];
                if (ctlg[v] <= 0) {
                    for (Index w : flowers[v]) {
                        if (ctlg[w] > 0) {
                            hit++;
                            tick++;
                            ctlg[v] = ctlg[w];
                            orphanEdges[j] = orphanEdges[static_cast<size_t>(ctlg[w]) - 1];
                        }
                    }
                }
            }
        }
        if (tick < orphanCount) {
            std::fprintf(stderr, "Error: seem to have missed some orphans\n");
        }
    }

    layoutVerts = intVerts;
    rimVerts = bdryList;
    v2indx.clear();
    indx2v.clear();

    return nodeCount;
}

Index Packer::pruneComplex() {
    if (orphanCount == 0) {
        return 0;
    }

    // Renamed from pruneComplex.m's local 'v2indx'/'indx2v' to
    // 'oldToNew'/'newToOld': an unqualified 'v2indx' inside this member
    // function would otherwise mean the Packer::v2indx *member*, which
    // serves an unrelated purpose (sparse-matrix layout indexing, set up by
    // indxMatrices()) -- MATLAB's obj.v2indx vs. a same-named local are
    // different namespaces, but C++ has no such distinction.
    const Index newBdryCount = static_cast<Index>(bdryList.size()) - 1;
    const Index newNodeCount = intCount + newBdryCount;
    std::vector<Index> oldToNew(static_cast<size_t>(nodeCount) + 1, 0);
    std::vector<Index> newToOld(static_cast<size_t>(newNodeCount) + 1, 0);
    for (Index j = 1; j <= intCount; ++j) {
        Index v = intVerts[static_cast<size_t>(j) - 1];
        newToOld[j] = v;
        oldToNew[v] = j;
    }
    for (Index j = 1; j <= newBdryCount; ++j) {
        Index w = bdryList[static_cast<size_t>(j) - 1];
        newToOld[intCount + j] = w;
        oldToNew[w] = intCount + j;
    }

    // ---- fix combinatorics ----
    std::vector<Index> newVNum(static_cast<size_t>(newNodeCount) + 1, 0);
    std::vector<Scalar> newVAims(static_cast<size_t>(newNodeCount) + 1, 0.0);
    std::vector<std::vector<Index>> newFlowers(static_cast<size_t>(newNodeCount) + 1);

    // Interior vertices: an interior vertex can never neighbor an orphan (by
    // definition of orphan -- cut off from the interior component), so every
    // petal survives; just renumber in place.
    for (Index nv = 1; nv <= intCount; ++nv) {
        Index v = intVerts[static_cast<size_t>(nv) - 1];
        newVNum[nv] = vNum[v];
        newVAims[nv] = vAims[v];
        const auto& flower = flowers[v];
        std::vector<Index> newflower(flower.size());
        for (size_t j = 0; j < flower.size(); ++j) {
            newflower[j] = oldToNew[flower[j]];
        }
        newFlowers[nv] = std::move(newflower);
    }

    // Boundary vertices: an orphan-adjacent petal is dropped; the flower is
    // rotated to start at its first surviving (non-orphan) petal so the
    // trimmed run stays contiguous, matching pruneComplex.m's own approach.
    for (Index nv = intCount + 1; nv <= newNodeCount; ++nv) {
        Index v = newToOld[nv];
        Index num = vNum[v];
        const auto& flower = flowers[v]; // size num+1 (open: front != back)

        Index tick = 0;
        Index spot = 0; // 1-indexed position of first surviving petal; 0 = none found yet
        while (spot == 0 && tick < num) {
            Index w = flower[static_cast<size_t>(tick)];
            if (oldToNew[w] > 0) {
                spot = tick + 1;
            } else {
                tick++;
            }
        }
        // Every vertex in bdryList is, by construction (complexCount()'s BFS
        // puts it there precisely because it neighbors some interior
        // vertex), guaranteed to have at least one surviving petal -- so
        // spot==0 here would mean an actual internal-consistency bug rather
        // than data pruneComplex.m was ever meant to handle (a literal port
        // would index flower(0), undefined in C++ and an error in MATLAB
        // too). Fail loudly instead of silently producing a corrupt Packer.
        if (spot == 0) {
            std::fprintf(stderr,
                         "Packer::pruneComplex: boundary vertex %d has no surviving "
                         "(non-orphan) petal; internal inconsistency\n", v);
            return 0;
        }

        std::vector<Index> newflower;
        newflower.reserve(static_cast<size_t>(num + 2 - spot));
        for (Index j = spot; j <= num + 1; ++j) {
            Index w = flower[static_cast<size_t>(j) - 1];
            Index nw = oldToNew[w];
            if (nw != 0) {
                newflower.push_back(nw);
            } else {
                break;
            }
        }
        newVNum[nv] = static_cast<Index>(newflower.size()) - 1;
        newVAims[nv] = -1.0;
        newFlowers[nv] = std::move(newflower);
    }

    // ---- organize combinatorial stuff ----
    const Index cutCount = nodeCount - newNodeCount;
    nodeCount = newNodeCount;
    orphanVerts.clear();
    orphanCount = 0;
    flowers = std::move(newFlowers);
    vAims = std::move(newVAims);
    vNum = std::move(newVNum); // harmless bookkeeping so vNum isn't
                                // momentarily stale; complexCount() below
                                // recomputes it from 'flowers' regardless,
                                // same as it always does.

    // translate int and bdry lists ('bdryList' is already in cclw order)
    for (Index j = 0; j < intCount; ++j) {
        intVerts[static_cast<size_t>(j)] = oldToNew[intVerts[static_cast<size_t>(j)]];
    }
    for (size_t j = 0; j < bdryList.size(); ++j) {
        bdryList[j] = oldToNew[bdryList[j]];
    }
    alpha = oldToNew[alpha];
    gamma = bdryList.front();

    // main organization
    //
    // Bug fix vs. pruneComplex.m: the source never checks complex_count()'s
    // return value here either -- see complexCount()'s own updated doc
    // comment for why silently continuing past a malformed complex is
    // worse than a clean, catchable error.
    if (complexCount() < 0) {
        throw std::runtime_error(
            "Packer::pruneComplex: complexCount() failed after pruning -- the "
            "trimmed complex has a malformed boundary (see the preceding stderr "
            "diagnostic)");
    }

    // ---- translate 'vlist' ----
    // pruneComplex.m reads 'obj.vlist(j)' here, AFTER 'obj.vlist' was just
    // cleared two lines earlier in the source -- real undefined behavior
    // (out-of-bounds access) in a literal C++ port, not just a MATLAB
    // quirk to preserve. Fixed to read from a saved local copy instead.
    std::vector<Index> savedVlist = vlist;
    vlist.clear();
    for (Index v : savedVlist) {
        Index nv = oldToNew[v];
        if (nv > 0) {
            vlist.push_back(nv);
        }
    }

    // ---- translate centers/radii ----
    std::vector<Scalar> newOrigRadii(static_cast<size_t>(newNodeCount) + 1, 0.0);
    std::vector<Complex> newOrigCenters(static_cast<size_t>(newNodeCount) + 1, Complex(0.0, 0.0));
    std::vector<Scalar> newRadii(static_cast<size_t>(newNodeCount) + 1, 0.0);
    std::vector<Complex> newCenters(static_cast<size_t>(newNodeCount) + 1, Complex(0.0, 0.0));
    std::vector<Scalar> newLocalRadii(static_cast<size_t>(newNodeCount) + 1, 0.0);
    std::vector<Complex> newLocalCenters(static_cast<size_t>(newNodeCount) + 1, Complex(0.0, 0.0));
    const bool haveOrig = !origRadii.empty() && !origCenters.empty();
    for (Index nv = 1; nv <= newNodeCount; ++nv) {
        Index v = newToOld[nv];
        if (haveOrig) {
            newOrigRadii[nv] = origRadii[v];
            // pruneComplex.m has 'origCenters(nv)=obj.origCenters(nv)' here,
            // almost certainly a typo for 'obj.origCenters(v)' -- it should
            // read from the OLD numbering, exactly like origRadii on the
            // line right above it and every other line in this loop. Ported
            // as the evident intent (a literal port would silently read
            // whichever old-numbered vertex happens to share nv's *new*
            // index -- wrong data, not just a crash, which is exactly the
            // kind of bug this project's port fixes rather than preserves).
            newOrigCenters[nv] = origCenters[v];
        }
        newRadii[nv] = radii[v];
        newCenters[nv] = centers[v];
        newLocalRadii[nv] = localradii[v];
        newLocalCenters[nv] = localcenters[v];
    }
    origRadii = std::move(newOrigRadii);
    origCenters = std::move(newOrigCenters);
    radii = std::move(newRadii);
    centers = std::move(newCenters);
    localradii = std::move(newLocalRadii);
    localcenters = std::move(newLocalCenters);

    // ---- outdated stuff ----
    angsumMonitor.clear();
    l2Monitor.clear();
    visErrMonitor.clear();
    ticMonitor.clear();
    corners.clear();
    sides.clear();

    std::fprintf(stdout, "Packing %s was pruned of %d orphan vertices\n", fileName.c_str(),
                 cutCount);

    return cutCount;
}

void Packer::indxMatrices(const std::vector<Index>& varlist) {
    if (varlist.empty()) {
        if (layoutVerts.empty()) {
            layoutVerts = intVerts;
            rimVerts = bdryList;
        }
    } else {
        std::vector<int> status(static_cast<size_t>(nodeCount) + 1, 0);
        if (hes == Geometry::Spherical) {
            for (Index w : bdryList) status[w] = 2;
        }

        layoutVerts.clear();
        for (Index v : varlist) {
            const auto& flower = flowers[v];
            if (flower.front() == flower.back() && status[v] <= 0) {
                for (Index w : flower) {
                    if (status[w] == 0) status[w] = -1;
                }
                layoutVerts.push_back(v);
                status[v] = 1;
            }
        }

        rimVerts.clear();
        for (Index v = 1; v <= nodeCount; ++v) {
            if (status[v] < 0) rimVerts.push_back(v);
        }

        if (static_cast<Index>(layoutVerts.size()) == intCount) {
            layoutVerts = intVerts;
            rimVerts = bdryList;
        }
    }

    const Index lolong = static_cast<Index>(layoutVerts.size());

    v2indx.assign(static_cast<size_t>(nodeCount) + 1, 0);
    indx2v.assign(1, 0); // indx2v[0] unused; entries pushed from index 1
    for (Index j = 0; j < lolong; ++j) {
        indx2v.push_back(layoutVerts[j]);
        v2indx[layoutVerts[j]] = static_cast<Index>(indx2v.size()) - 1;
    }
    for (size_t j = 0; j + 1 < rimVerts.size(); ++j) { // rimVerts(1:end-1)
        Index w = rimVerts[j];
        indx2v.push_back(w);
        v2indx[w] = static_cast<Index>(indx2v.size()) - 1;
    }

    if (varlist.empty()) {
        Index minv2indx = nodeCount > 0 ? v2indx[1] : 0;
        for (Index v = 1; v <= nodeCount; ++v) minv2indx = std::min(minv2indx, v2indx[v]);
        if (minv2indx == 0) {
            std::fprintf(stderr, "Warning: not all vertices have interior neighbors;\n");
            std::fprintf(stderr, "  there may be packing or layout problems.\n");
        }
    }

    Index ijCount = 0;
    for (Index k = 1; k <= lolong; ++k) {
        ijCount += vNum[indx2v[k]] + 1;
    }

    tranI.assign(static_cast<size_t>(ijCount), 0);
    tranJ.assign(static_cast<size_t>(ijCount), 0);
    tranJindx.assign(static_cast<size_t>(ijCount), 0);
    rhsI.assign(static_cast<size_t>(ijCount), 0);
    rhsJ.assign(static_cast<size_t>(ijCount), 0);
    rhsJindx.assign(static_cast<size_t>(ijCount), 0);

    Index kj = 0; // 0-based running count of tran entries
    Index kw = 0; // 0-based running count of rhs entries
    for (Index k = 1; k <= lolong; ++k) {
        Index v = indx2v[k];

        // diagonal entry first
        tranI[kj] = k;
        tranJ[kj] = k;
        tranJindx[kj] = -1; // edge to self
        kj++;

        const auto& flower = flowers[v];
        Index num = static_cast<Index>(flower.size()) - 1;
        for (Index j = 1; j <= num; ++j) {
            Index w = flower[j - 1]; // flower(j) in MATLAB, 1-indexed -> 0-indexed vector
            // Bug fix / hardening vs. indxMatrices.m: 'v2indx(w)<=lolong' is
            // the literal source condition, but v2indx[w]==0 means w was
            // never assigned ANY index -- it's neither a layout vertex nor a
            // rim vertex. That should be impossible for a well-formed
            // complex (every neighbor of an interior vertex is either
            // interior itself or an immediate boundary neighbor, and both
            // categories always get indexed above), but 0 still satisfies
            // '<= lolong' for any lolong>=0, so a literal port silently
            // mis-files it as an interior neighbor with column index 0.
            // That corrupts the sparse-triplet lists layoutCenters() builds
            // from tranJ/rhsJ: a triplet column of 0 there becomes -1 after
            // the 1-indexed-to-0-indexed shift, which Eigen's
            // setFromTriplets() has no obligation to bounds-check -- in a
            // Release build that reads as heap corruption / an
            // unpredictable access violation somewhere downstream in
            // riffle(), not a clean, diagnosable error at the actual point
            // of corruption. Fail loudly and immediately instead, with
            // enough detail to identify which vertex/complex is malformed
            // (most likely cause: complexCount() produced an incomplete
            // bdryList -- see that function's own updated doc comment).
            if (v2indx[w] == 0) {
                throw std::runtime_error(
                    "Packer::indxMatrices: vertex " + std::to_string(w) +
                    " (a neighbor of interior vertex " + std::to_string(v) +
                    ") was never assigned a layout or rim index -- the complex's "
                    "boundary/interior bookkeeping is inconsistent (see complexCount())");
            }
            if (v2indx[w] <= lolong) {
                tranI[kj] = k;
                tranJ[kj] = v2indx[w];
                tranJindx[kj] = j;
                kj++;
            } else {
                rhsI[kw] = k;
                rhsJ[kw] = v2indx[w] - lolong;
                rhsJindx[kw] = j;
                kw++;
            }
        }
    }

    tranIJcount = kj;
    rhsIJcount = kw;
    tranI.resize(static_cast<size_t>(kj));
    tranJ.resize(static_cast<size_t>(kj));
    tranJindx.resize(static_cast<size_t>(kj));
    rhsI.resize(static_cast<size_t>(kw));
    rhsJ.resize(static_cast<size_t>(kw));
    rhsJindx.resize(static_cast<size_t>(kw));
}

} // namespace gopack
