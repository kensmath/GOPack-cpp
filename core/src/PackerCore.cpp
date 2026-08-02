#include "gopack/Packer.h"

#include <algorithm>
#include <cstdio>
#include <deque>
#include <random>

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
        if (nextb != firstbdry || bdryList.back() != bdryList.front()) {
            std::fprintf(stderr, "Error forming bdryList\n");
        }
        if (bdryList.size() >= maxTick) {
            std::fprintf(stderr, "Error in bdryList, too long\n");
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
    throw NotImplementedError(
        "Packer::pruneComplex: orphan-vertex removal has not been ported yet "
        "(pruneComplex.m). This is only exercised by the random rectangle/"
        "square generators, not by loading real packing files, so it was "
        "deferred out of the core max-pack port.");
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
