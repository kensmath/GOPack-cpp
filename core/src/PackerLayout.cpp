#include "gopack/Packer.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>

namespace gopack {

using geom::cosAngle;
using geom::cosCorner;

int Packer::setMode(int mdIn, const std::vector<Index>& crns, const std::vector<Scalar>& angs) {
    if (nodeCount < 1) {
        std::fprintf(stderr, "GOpacker seems to have no packing.\n");
        mode = -1;
        return -1;
    }

    int m = std::max(mdIn, 1);
    if (m < 1 || (m > 2 && m != 4)) {
        std::fprintf(stderr,
                     "Mode choices: 1 = max packing, 2 = polygonal packing, "
                     "4 = orthopack (euclidean disc packing, boundary orthogonal to unit circle).\n");
        return -1;
    }
    int md = m;

    if (m > 1 && hes == Geometry::Spherical) {
        std::fprintf(stderr, "Mode must be %s for a sphere.\n", kPackModes[0]);
        md = 1;
        m = 1;
    }

    if (m == 1) {
        vAims.assign(static_cast<size_t>(nodeCount) + 1, kTwoPi);
        for (Index w = 1; w <= nodeCount; ++w) {
            if (bdryFlags[w] != 0) vAims[w] = -1.0;
        }
        mode = md;
        std::fprintf(stderr, "Mode is set to %s\n", kPackModes[mode - 1]);
        return md;
    }

    // ---- m == 4: orthopack (new; see setOrthoCenters()'s doc comment and
    // the "MODE 4" module comment in Packer.h) ----
    if (m == 4) {
        // "Triangulates a disc" check (see setMode()'s own doc comment):
        // complexCount() has already rejected a malformed/disconnected
        // boundary walk by the time any caller could reach setMode(), so the
        // one thing left to check here is that there IS a boundary at all --
        // a closed/spherical complex (bdryCount==0) has no boundary circles
        // to lay out orthogonally in the first place.
        if (bdryCount <= 0) {
            std::fprintf(stderr,
                         "orthopack requires the complex to triangulate a disc (a nonempty "
                         "boundary) -- this complex is closed (spherical), with no boundary.\n");
            mode = -1;
            return -1;
        }

        // Deliberate deviation from GOPacker.m, same rationale as mode 2
        // below: orthopack always produces a euclidean packing (boundary
        // circles orthogonal to the unit circle is inherently a euclidean
        // notion), regardless of what geometry the input was read as.
        hes = Geometry::Euclidean;

        // Same vAims convention as mode 1 (interior target angle sum 2*pi;
        // boundary vertices get the -1 "free" sentinel setEffective() already
        // branches on) -- orthopack's boundary radii, like mode 1's, are
        // determined by the packing process (Steps A/B/C) rather than by a
        // fixed target, so they need the same treatment mode 1 gets, not
        // mode 2's fixed corner/side angle targets.
        vAims.assign(static_cast<size_t>(nodeCount) + 1, kTwoPi);
        for (Index w = 1; w <= nodeCount; ++w) {
            if (bdryFlags[w] != 0) vAims[w] = -1.0;
        }
        mode = md;
        std::fprintf(stderr, "Mode is set to %s\n", kPackModes[mode - 1]);
        return md;
    }

    // ---- m == 2: polygonal packing ----

    // Deliberate deviation from GOPacker.m: force euclidean geometry for
    // polygonal/rectangle mode, even if the packing was originally read as
    // hyperbolic or spherical. GOPack always computes internally in
    // euclidean coordinates regardless of 'hes' (see the module comment in
    // Geometry.h) -- 'hes' only controls whether readpack()/writepack()
    // convert to/from hyperbolic or spherical circle data at the file
    // boundary. A polygon/rectangle boundary is inherently a euclidean
    // shape (straight sides, corner angles measured in the plane), so
    // re-interpreting its output as hyperbolic or spherical on write
    // wouldn't produce a meaningful hyperbolic/spherical polygon -- it
    // would just silently warp already-euclidean coordinates through a
    // conversion that was never intended for them. Setting hes here means
    // writepack() leaves the (already-euclidean) radii/centers alone
    // instead of applying that conversion.
    hes = Geometry::Euclidean;

    // GOPacker.m's `obj.vAims(v)=val` assignments below auto-grow MATLAB's
    // array if it's currently shorter than `v` -- std::vector::operator[]
    // has no such safety net, so writing to an under-sized vAims here is
    // undefined behavior (in practice: silent heap corruption in a Release
    // build, not a clean crash). readpack() always leaves vAims sized
    // nodeCount+1 before a caller could reach mode 2, but guarantee it here
    // too so setMode(2, ...) is safe to call on its own, matching what the
    // m==1 branch above already does unconditionally.
    if (vAims.size() < static_cast<size_t>(nodeCount) + 1) {
        vAims.resize(static_cast<size_t>(nodeCount) + 1, kTwoPi);
    }

    Index sideN = 4; // default number of sides for a polygon
    std::vector<Index> cornersLocal; // local listing of corners, if given

    // GOPacker.m distinguishes "crns wasn't passed at all" (nargin<3) from
    // "crns was passed as an explicit empty array" (nargin>=3, cln==0) --
    // C++ has no such distinction for a by-value vector default argument,
    // so an empty 'crns' here is treated as "not given", falling straight
    // through to the vlist/random corner-selection below, rather than
    // MATLAB's stricter immediate "empty or too short" error. This is the
    // more useful behavior for a caller (setMode(2) alone should mean "you
    // figure out the corners"), and it only changes behavior for an
    // explicitly-empty-but-passed argument, which C++ can't distinguish
    // from "omitted" anyway.
    if (!crns.empty()) {
        Index cln = static_cast<Index>(crns.size());
        if (cln == 2) {
            std::fprintf(stderr, "Error: list of corners is empty or too short\n");
            mode = -1;
            return -1;
        }
        if (cln == 1) {
            if (crns[0] > 2) {
                sideN = crns[0];
                if (!vlist.empty()) {
                    if (static_cast<Index>(vlist.size()) > sideN) {
                        vlist.resize(static_cast<size_t>(sideN)); // truncate
                    } else if (static_cast<Index>(vlist.size()) < sideN) {
                        vlist.clear(); // discard vlist
                    }
                }
            }
        } else {
            sideN = cln;
            cornersLocal.assign(static_cast<size_t>(sideN), 0);
            for (Index k = 0; k < cln; ++k) {
                Index cv = crns[static_cast<size_t>(k)];
                if (bdryFlags[cv] != 1) {
                    std::fprintf(stderr,
                                 "Error: Your given \"corner\" %d is not a boundary vertex\n", cv);
                    mode = -1;
                    return -1;
                }
                cornersLocal[static_cast<size_t>(k)] = cv;
            }
        }

        if (!angs.empty()) {
            if (static_cast<Index>(angs.size()) != sideN) {
                std::fprintf(stderr, "Numbers of corner vertices and angles do not match");
                mode = -1;
                return -1;
            }
            Scalar sum = 0.0;
            for (Scalar a : angs) sum += a;
            if (std::abs(sideN * kPi - sum - kTwoPi) > 0.01) {
                std::fprintf(stderr,
                             "Corner angles are not consistent with polygon turning angles");
                mode = -1;
                return -1;
            }
        }
    }

    // Corners not given? take bdry vertices in 'vlist' and infer sideN.
    if (cornersLocal.empty() && vlist.size() >= 3) {
        Index vln = static_cast<Index>(vlist.size());
        cornersLocal.assign(static_cast<size_t>(vln), 0);
        Index mctn = 0;
        for (Index k = 0; k < vln; ++k) {
            if (bdryFlags[vlist[static_cast<size_t>(k)]] != 0) {
                cornersLocal[static_cast<size_t>(mctn)] = vlist[static_cast<size_t>(k)];
                mctn++;
            }
        }
        if (mctn < 3) { // didn't get enough? default to random below
            cornersLocal.clear();
        } else {
            sideN = mctn;
            cornersLocal.resize(static_cast<size_t>(sideN)); // trim to right length
        }
    }

    // Still no corners? choose 'sideN' corners pseudo-randomly around bdryList.
    if (cornersLocal.empty()) {
        Index bl = static_cast<Index>(bdryList.size()); // closed list: bdryCount+1
        if (bl < 3) {
            std::fprintf(stderr, "The boundary has only %d vertices\n", bl);
            mode = -1;
            return -1;
        } else if (bl == 3) { // triangle
            sideN = 3;
        } else if (bl < sideN) { // bl is maximum number of sides
            sideN = bl;
        }

        cornersLocal.assign(static_cast<size_t>(sideN), 0);
        Index fth = bl / sideN;
        static thread_local std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<Index> dist(1, bl);
        Index sd = dist(rng);
        for (Index e = 0; e < sideN; ++e) {
            cornersLocal[static_cast<size_t>(e)] = bdryList[static_cast<size_t>(sd) - 1];
            sd = 1 + (sd + fth) % bl;
        }
        vlist = cornersLocal;
        std::fprintf(stderr, "Corners not provided, so %d were chosen randomly\n", sideN);
    }

    // check length
    sideN = static_cast<Index>(cornersLocal.size());
    if (sideN < 3) {
        std::fprintf(stderr, "Setting polygon mode requires at least 3 corners\n");
        mode = -1;
        return -1;
    }

    // set all bdry aims to pi, equal corner target angles as default
    for (Index i = 1; i <= bdryCount; ++i) {
        vAims[static_cast<size_t>(bdryList[static_cast<size_t>(i) - 1])] = kPi;
    }
    for (Index i = 0; i < sideN; ++i) {
        vAims[static_cast<size_t>(cornersLocal[static_cast<size_t>(i)])] =
            kPi * (1.0 - 2.0 / sideN);
    }

    if (!angs.empty()) { // 'angs' were specified
        for (Index i = 0; i < sideN; ++i) {
            vAims[static_cast<size_t>(cornersLocal[static_cast<size_t>(i)])] =
                angs[static_cast<size_t>(i)];
        }
    }

    // determine indices in bdryList so we can get cclw order
    std::vector<Index> bdryIndx(static_cast<size_t>(sideN), 0);
    for (Index i = 0; i < sideN; ++i) {
        Index cnr = cornersLocal[static_cast<size_t>(i)];
        for (Index j = 1; j <= bdryCount; ++j) {
            if (bdryList[static_cast<size_t>(j) - 1] == cnr) {
                bdryIndx[static_cast<size_t>(i)] = j;
                break;
            }
        }
        if (bdryIndx[static_cast<size_t>(i)] == 0) {
            std::fprintf(stderr, "Vert %d is not a bdry vertex\n", cnr);
            mode = -1;
            return -1;
        }
    }

    // The source's "put 'cornangs' in counterclockwise order" block
    // (tmpcorners=corners; corners=zeros(...); corners(j)=bdryList(bdryIndx(j)))
    // recomputes a local 'corners' variable that is never read again, and
    // does so using an unsorted bdryIndx (its `sort(bdryIndx);` call never
    // captures the result) -- so that whole block is dead code, faithfully
    // omitted here. The only thing it produces that's used later is
    // 'upright' = the first requested corner.
    Index upright = cornersLocal[0];

    std::vector<Index> sortedBdryIndx = bdryIndx;
    std::sort(sortedBdryIndx.begin(), sortedBdryIndx.end());

    Index offset = 1;
    for (Index i = 2; i <= sideN; ++i) {
        if (sortedBdryIndx[static_cast<size_t>(i) - 1] == upright) offset = i;
    }

    corners.assign(static_cast<size_t>(sideN), 0);
    std::vector<Index> crnrIndices(static_cast<size_t>(sideN), 0);
    for (Index i = 0; i < sideN; ++i) {
        Index k = 1 + (offset - 1 + i) % sideN;
        crnrIndices[static_cast<size_t>(i)] = sortedBdryIndx[static_cast<size_t>(k) - 1];
        corners[static_cast<size_t>(i)] =
            bdryList[static_cast<size_t>(crnrIndices[static_cast<size_t>(i)]) - 1];
    }

    // set up 'sides': for each corner, walk the closed bdryList forward
    // (by bdryList-position, wrapping mod bdryCount) until the next corner
    // is reached, collecting every boundary vertex along that side.
    sides.assign(static_cast<size_t>(sideN), {});
    for (Index i = 1; i <= sideN; ++i) {
        Index cornerindx = crnrIndices[static_cast<size_t>(i) - 1];
        Index j = 1 + (i % sideN);
        Index nextcorner = crnrIndices[static_cast<size_t>(j) - 1];

        std::vector<Index> sideIndices(static_cast<size_t>(bdryCount), 0);
        sideIndices[0] = crnrIndices[static_cast<size_t>(i) - 1];
        Index tick = 1;
        while (sideIndices[static_cast<size_t>(tick) - 1] != nextcorner) {
            Index k = (cornerindx + tick - 1) % bdryCount + 1;
            tick++;
            if (tick > bdryCount) {
                std::fprintf(stderr, "error in getting sides\n");
                mode = -1;
                return -1;
            }
            sideIndices[static_cast<size_t>(tick) - 1] = k;
        }

        std::vector<Index> side(static_cast<size_t>(tick), 0);
        for (Index jj = 1; jj <= tick; ++jj) {
            side[static_cast<size_t>(jj) - 1] =
                bdryList[static_cast<size_t>(sideIndices[static_cast<size_t>(jj) - 1]) - 1];
        }
        sides[static_cast<size_t>(i) - 1] = std::move(side);
    }

    mode = m;
    md = m;
    std::fprintf(stdout, "Mode is \"%s\", corner vertices are:", kPackModes[mode - 1]);
    for (Index c : corners) std::fprintf(stdout, " %d", c);
    std::fprintf(stdout, "\n");
    return md;
}

void Packer::layoutBdry() {
    if (mode == 2) {
        setPolyCenters();
        return;
    }
    if (mode == 4) {
        setOrthoCenters();
        return;
    }
    setHoroCenters();
}

void Packer::setRectCenters() {
    // Assumes mode==2, corners/sides already set by setMode(2, ...), and
    // exactly 4 sides (setPolyCenters only calls this when num_sides==4
    // and all four corner aims are within 1e-5 of a right angle).
    std::vector<Scalar> sidelengths(4, 0.0);
    for (int i = 0; i < 4; ++i) {
        const auto& side = sides[static_cast<size_t>(i)];
        Index n = static_cast<Index>(side.size());
        Scalar longv = localradii[static_cast<size_t>(side[0])];
        longv += localradii[static_cast<size_t>(side[static_cast<size_t>(n) - 1])];
        for (Index j = 2; j <= n - 1; ++j) {
            longv += 2.0 * localradii[static_cast<size_t>(side[static_cast<size_t>(j) - 1])];
        }
        sidelengths[static_cast<size_t>(i)] = longv;
    }

    Scalar width = (sidelengths[0] + sidelengths[2]) / 2.0;
    Scalar height = (sidelengths[1] + sidelengths[3]) / 2.0;

    Scalar aspect = height / width;
    Scalar factor = 2.0 * (aspect + 1.0) / (width + height);
    for (Index v = 1; v <= nodeCount; ++v) localradii[static_cast<size_t>(v)] *= factor;
    for (auto& s : sidelengths) s *= factor;

    // NOTE ("?????" in the source): despite the source comment describing
    // "lowerleft (-aspect,-1), upper right (aspect,1)", the actual corner
    // points below place the rectangle's real part in {+1,-1} and its
    // imaginary part in {+aspect,-aspect} -- ported literally as written,
    // not as commented.
    const Complex crnpt[4] = {Complex(1.0, aspect), Complex(-1.0, aspect),
                               Complex(-1.0, -aspect), Complex(1.0, -aspect)};
    const Complex edgedir[4] = {Complex(-1.0, 0.0), Complex(0.0, -1.0), Complex(1.0, 0.0),
                                 Complex(0.0, 1.0)};
    const Scalar slength[4] = {2.0, 2.0 * aspect, 2.0, 2.0 * aspect};

    for (int k = 0; k < 4; ++k) {
        Scalar sidefactor = slength[static_cast<size_t>(k)] / sidelengths[static_cast<size_t>(k)];
        const auto& side = sides[static_cast<size_t>(k)];
        Index n = static_cast<Index>(side.size());
        Scalar prev = localradii[static_cast<size_t>(corners[static_cast<size_t>(k)])];
        Complex spot = crnpt[static_cast<size_t>(k)];
        localcenters[static_cast<size_t>(side[0])] = spot;
        // Only walks to n-2 (not n-1): the side's last vertex is the next
        // corner, whose exact position is set directly as crnpt[k+1] at
        // the top of the next iteration rather than accumulated here.
        for (Index i = 1; i <= n - 2; ++i) {
            Scalar next = localradii[static_cast<size_t>(side[static_cast<size_t>(i)])];
            spot += sidefactor * edgedir[static_cast<size_t>(k)] * (prev + next);
            localcenters[static_cast<size_t>(side[static_cast<size_t>(i)])] = spot;
            prev = next;
        }
    }
}

void Packer::setPolyCenters() {
    if (mode != 2) {
        std::fprintf(stderr, "setPolyCenters: mode must be %s\n", kPackModes[1]);
        return;
    }

    Index numSides = static_cast<Index>(corners.size());
    if (numSides < 3 || numSides != static_cast<Index>(sides.size()) || vAims.empty()) {
        // Faithful to the source: setPolyCenters.m prints this diagnostic
        // but has no early return here, so execution falls through even
        // when the precondition fails. In practice this is unreachable
        // when corners/sides came from setMode(2, ...), which guarantees
        // numSides>=3 and a matching 'sides' array.
        std::fprintf(stderr, "setPolyCenters: must have corners, sides, and aims.\n");
    }

    // check for rectangle first
    if (numSides == 4) {
        Scalar benderror = 0.0;
        for (Index j = 0; j < 4; ++j) {
            benderror +=
                std::abs(vAims[static_cast<size_t>(corners[static_cast<size_t>(j)])] - kPi / 2.0);
        }
        if (benderror <= 0.00001) { // standard situation, right angles
            setRectCenters();
            return;
        }
    } // else, handled below by the n=even routines

    // compute side lengths using 'localradii'
    std::vector<Scalar> sidelengths(static_cast<size_t>(numSides), 0.0);
    Scalar fullLength = 0.0;
    std::vector<Scalar> targetLength(static_cast<size_t>(numSides), 1.0);
    for (Index i = 0; i < numSides; ++i) {
        const auto& side = sides[static_cast<size_t>(i)];
        Index n = static_cast<Index>(side.size());
        Scalar longv = localradii[static_cast<size_t>(side[0])];
        longv += localradii[static_cast<size_t>(side[static_cast<size_t>(n) - 1])];
        for (Index j = 2; j <= n - 1; ++j) {
            longv += 2.0 * localradii[static_cast<size_t>(side[static_cast<size_t>(j) - 1])];
        }
        sidelengths[static_cast<size_t>(i)] = longv;
        fullLength += longv;
    }
    Index halfn = numSides / 2;

    if (numSides == 3) {
        // triangle: solve triangle to set 'targetLength's
        Scalar opp1 = vAims[static_cast<size_t>(corners[2])];
        Scalar opp2 = vAims[static_cast<size_t>(corners[0])];
        if (opp1 <= 0.0 || opp2 <= 0.0 || (opp1 + opp2) >= kPi) {
            std::fprintf(stderr, "setPolyCenters: error in triangles angle aims\n");
            return;
        }
        Scalar opp3 = kPi - (opp1 + opp2); // ensure angles sum to pi
        // law of sines gives desired proportions of side lengths
        targetLength[0] = 1.0;
        targetLength[1] = std::sin(opp2) / std::sin(opp1);
        targetLength[2] = std::sin(opp3) / std::sin(opp1);
        Scalar lensum = targetLength[0] + targetLength[1] + targetLength[2];
        Scalar factor = lensum / fullLength;
        for (auto& r : localradii) r *= factor;
        for (auto& s : sidelengths) s *= factor;
    } else if (halfn * 2 == numSides) {
        // polygon, n even: pair up opposite sides, target total length ~2*pi
        Scalar factor = 6.0 / fullLength;
        for (auto& r : localradii) r *= factor;
        for (auto& s : sidelengths) s *= factor;
        for (Index j = 1; j <= halfn; ++j) {
            targetLength[static_cast<size_t>(j) - 1] =
                (sidelengths[static_cast<size_t>(j) - 1] +
                 sidelengths[static_cast<size_t>(halfn + j) - 1]) /
                2.0;
            targetLength[static_cast<size_t>(halfn + j) - 1] =
                targetLength[static_cast<size_t>(j) - 1];
        }
    } else {
        // polygon, n odd: sides target length 2*sin(pi/n) (regular n-gon in unit disc)
        Scalar spn = 2.0 * std::sin(kPi / numSides);
        Scalar factor = numSides * spn / fullLength;
        for (auto& r : localradii) r *= factor;
        for (auto& s : sidelengths) s *= factor;
        for (Index j = 0; j < numSides; ++j) targetLength[static_cast<size_t>(j)] = spn;
    }

    // ------ Step 1: lay out using 'targetLength's.
    //   num_sides odd: first corner at z=i, bisected by imaginary axis,
    //     edge down to left.
    //   num_sides even: first corner at 1+i, first edge horizontal to left.
    std::vector<Scalar> edgeArg(static_cast<size_t>(numSides), 1.0);
    edgeArg[0] = kPi;
    if (halfn * 2 != numSides) { // odd?
        edgeArg[0] = kPi + (kPi - vAims[static_cast<size_t>(corners[0])]) / 2.0;
    }
    for (Index j = 2; j <= numSides; ++j) {
        edgeArg[static_cast<size_t>(j) - 1] =
            edgeArg[static_cast<size_t>(j) - 2] + kPi -
            vAims[static_cast<size_t>(corners[static_cast<size_t>(j) - 1])];
    }
    std::vector<Complex> edgedir(static_cast<size_t>(numSides), Complex(1.0, 0.0));
    for (Index j = 0; j < numSides; ++j) {
        edgedir[static_cast<size_t>(j)] = std::exp(Complex(0.0, edgeArg[static_cast<size_t>(j)]));
    }

    // Set first corner, then layout edges in turn; adjust sizes based on 'targetLength's.
    localcenters[static_cast<size_t>(corners[0])] = Complex(0.0, 1.0); // at z=i
    if (halfn * 2 == numSides) {                                      // even?
        localcenters[static_cast<size_t>(corners[0])] = Complex(1.0, 1.0); // at z=1+i
    }
    for (Index k = 1; k <= numSides; ++k) {
        Scalar sidefactor =
            targetLength[static_cast<size_t>(k) - 1] / sidelengths[static_cast<size_t>(k) - 1];
        const auto& side = sides[static_cast<size_t>(k) - 1];
        Index n = static_cast<Index>(side.size());
        Scalar prev = localradii[static_cast<size_t>(corners[static_cast<size_t>(k) - 1])];
        Complex spot = localcenters[static_cast<size_t>(corners[static_cast<size_t>(k) - 1])];
        localcenters[static_cast<size_t>(side[0])] = spot;
        for (Index i = 1; i <= n - 1; ++i) {
            Scalar next = localradii[static_cast<size_t>(side[static_cast<size_t>(i)])];
            spot += sidefactor * edgedir[static_cast<size_t>(k) - 1] * (prev + next);
            localcenters[static_cast<size_t>(side[static_cast<size_t>(i)])] = spot;
            prev = next;
        }
    }

    // Step 2: put average of corners at the origin.
    Complex centAvg(0.0, 0.0);
    for (Index j = 0; j < numSides; ++j) {
        centAvg += localcenters[static_cast<size_t>(corners[static_cast<size_t>(j)])];
    }
    centAvg /= static_cast<Scalar>(numSides);
    for (Index j = 1; j <= nodeCount; ++j) localcenters[static_cast<size_t>(j)] -= centAvg;

    // Step 3: scale.
    Scalar scalefactor;
    if (halfn * 2 == numSides) { // even?
        scalefactor = localcenters[static_cast<size_t>(corners[0])].real();
    } else {
        scalefactor = localcenters[static_cast<size_t>(corners[0])].imag();
    }
    for (auto& c : localcenters) c /= scalefactor;
    for (auto& r : localradii) r /= scalefactor;
}

Scalar Packer::getAspect() {
    if (mode != 2 || corners.size() != 4) {
        std::fprintf(stderr, "Aspect usage: should be mode 2 and have 4 corners\n");
        return -1.0;
    }
    Scalar top = std::abs(localcenters[static_cast<size_t>(corners[1])] -
                          localcenters[static_cast<size_t>(corners[0])]);
    Scalar rend = std::abs(localcenters[static_cast<size_t>(corners[3])] -
                           localcenters[static_cast<size_t>(corners[0])]);
    Scalar lend = std::abs(localcenters[static_cast<size_t>(corners[2])] -
                           localcenters[static_cast<size_t>(corners[1])]);
    Scalar bot = std::abs(localcenters[static_cast<size_t>(corners[3])] -
                          localcenters[static_cast<size_t>(corners[2])]);
    return (top + bot) / (rend + lend);
}

void Packer::setHoroCenters() {
    if (bdryCount <= 3) {
        Scalar s3 = std::sqrt(3.0);
        Scalar brad = s3 / (2 + s3);
        for (Index i = 0; i < 3 && i < bdryCount; ++i) {
            Index bv = bdryList[static_cast<size_t>(i)];
            localradii[bv] = brad;
            vAims[bv] = 0.0;
        }
        if (bdryCount >= 1) localcenters[bdryList[0]] = Complex(0.0, 1.0 - brad);
        if (bdryCount >= 2)
            localcenters[bdryList[1]] = (1 - brad) * Complex(-std::sqrt(3.0) / 2.0, -0.5);
        if (bdryCount >= 3)
            localcenters[bdryList[2]] = (1 - brad) * Complex(std::sqrt(3.0) / 2.0, -0.5);
        return;
    }

    // initial guess for R
    Scalar R = 0.0;
    Scalar minrad = 0.0;
    std::vector<Scalar> r(static_cast<size_t>(bdryCount) + 2, 0.0); // 1-indexed, closed
    for (Index j = 1; j <= bdryCount; ++j) {
        r[j] = localradii[bdryList[static_cast<size_t>(j) - 1]];
        if (r[j] > minrad) minrad = r[j];
        R += r[j];
    }
    r[bdryCount + 1] = r[1];
    R /= kPi;
    if (R < 2.0 * minrad) R = 3.0 * minrad;

    // Newton iteration to find R
    int trys = 0;
    bool keepon = true;
    while (keepon && trys < 100) {
        trys++;
        Scalar fvalue = -2.0 * kPi;
        Scalar fprime = 0.0;
        for (Index j = 1; j <= bdryCount; ++j) {
            Scalar Rrr = R - r[j] - r[j + 1];
            Scalar RRrr = R * Rrr;
            Scalar ab = r[j] * r[j + 1];
            fvalue += std::acos((RRrr - ab) / (RRrr + ab));
            fprime -= 1.0 * (R + Rrr) * std::sqrt(ab / RRrr) / (RRrr + ab);
        }

        Scalar newR = R - fvalue / fprime;
        if (newR < R / 2.0) newR = R / 2.0;
        if (newR > 2.0 * R) newR = 2.0 * R;
        if (std::abs(newR - R) < 0.00001) keepon = false;
        R = newR;
    }

    for (Index v = 1; v <= nodeCount; ++v) localradii[v] /= R;
    for (Index j = 1; j <= bdryCount + 1; ++j) r[j] /= R;

    Scalar r2 = r[1];
    localcenters[bdryList[0]] = Complex(0.0, 1.0 - r2);
    Scalar arg = kPi / 2.0;
    for (Index k = 2; k <= bdryCount; ++k) {
        Scalar r1 = r2;
        r2 = r[k];
        Scalar RRrr = 1.0 - r1 - r2;
        Scalar ab = r1 * r2;
        Scalar delta = std::acos((RRrr - ab) / (RRrr + ab));
        arg += delta;
        Scalar d = 1.0 - r2;
        localcenters[bdryList[static_cast<size_t>(k) - 1]] =
            Complex(d * std::cos(arg), d * std::sin(arg));
    }
}

// setOrthoCenters -- new, no MATLAB counterpart. See the doc comment on the
// declaration in Packer.h for the full derivation; this is the euclidean
// counterpart of setHoroCenters() immediately above (structured the same
// way -- gather boundary radii, Newton-solve for a common-circle radius R,
// normalize, walk bdryList placing centers), except boundary circles land
// orthogonal to the common (unit, after normalization) circle instead of
// internally tangent to it.
void Packer::setOrthoCenters() {
    if (bdryCount < 3) {
        std::fprintf(stderr, "setOrthoCenters: need at least 3 boundary vertices\n");
        return;
    }

    // gather boundary radii, closed (r[bdryCount+1] duplicates r[1], as in
    // setHoroCenters())
    Scalar sumR = 0.0;
    Scalar minrad = 0.0;
    std::vector<Scalar> r(static_cast<size_t>(bdryCount) + 2, 0.0); // 1-indexed, closed
    for (Index j = 1; j <= bdryCount; ++j) {
        r[j] = localradii[bdryList[static_cast<size_t>(j) - 1]];
        if (r[j] > minrad) minrad = r[j];
        sumR += r[j];
    }
    r[bdryCount + 1] = r[1];

    const Index n = bdryCount;
    const Scalar target = static_cast<Scalar>(n - 2) * kPi;

    // Initial guess for R: the root satisfies 0 < R < sumR/pi (see the
    // header doc comment), so start below that bound; minrad guards against
    // sumR/pi landing at (numerically) zero for degenerate all-but-zero
    // radii.
    Scalar R = 0.5 * sumR / kPi;
    if (!(R > 0.0) || R < 1e-9 * std::max(minrad, 1.0)) {
        R = (minrad > 0.0) ? 0.5 * minrad : 1.0;
    }

    // Newton iteration to solve sum_j 2*atan(R/r_j) = (n-2)*pi for R (see
    // the header doc comment for why this has a unique positive root).
    // Structured identically to setHoroCenters()'s own Newton loop just
    // above (same clamp-the-step/trys<100/convergence-tolerance shape), just
    // with this mode's own fvalue/fprime.
    int trys = 0;
    bool keepon = true;
    while (keepon && trys < 100) {
        trys++;
        Scalar fvalue = -target;
        Scalar fprime = 0.0;
        for (Index j = 1; j <= n; ++j) {
            fvalue += 2.0 * std::atan(R / r[j]);
            fprime += 2.0 * r[j] / (r[j] * r[j] + R * R);
        }
        Scalar newR = R - fvalue / fprime;
        if (newR < R / 2.0) newR = R / 2.0;
        if (newR > 2.0 * R) newR = 2.0 * R;
        if (std::abs(newR - R) < 1e-9 * std::max<Scalar>(R, 1.0)) keepon = false;
        R = newR;
    }

    // Normalize so the common orthogonal circle is exactly the unit circle:
    // scaling every radius (and R itself) by the same factor 1/R preserves
    // both the orthogonality identity d^2=R^2+r^2 and tangency (both are
    // similarity-invariant), so after this rescaling the common circle's
    // radius is exactly 1.
    for (Index v = 1; v <= nodeCount; ++v) localradii[v] /= R;
    for (Index j = 1; j <= bdryCount + 1; ++j) r[j] /= R;

    // d[j] = distance from the origin to boundary circle j's center, on the
    // now-unit-circle-normalized scale (R==1): d^2 = 1 + r_j^2.
    std::vector<Scalar> d(static_cast<size_t>(bdryCount) + 2, 0.0);
    for (Index j = 1; j <= bdryCount + 1; ++j) d[j] = std::sqrt(1.0 + r[j] * r[j]);

    // Walk cclw around bdryList, starting with bdryList[0] straight up
    // (matching setHoroCenters()'s own starting orientation), placing each
    // subsequent center at the cumulative angle from the law-of-cosines
    // formula above.
    localcenters[bdryList[0]] = Complex(0.0, d[1]);
    Scalar arg = kPi / 2.0;
    for (Index k = 2; k <= bdryCount; ++k) {
        Scalar r1 = r[static_cast<size_t>(k) - 1];
        Scalar r2 = r[static_cast<size_t>(k)];
        Scalar cosTheta = (1.0 - r1 * r2) / (d[static_cast<size_t>(k) - 1] * d[static_cast<size_t>(k)]);
        cosTheta = std::max(-1.0, std::min(1.0, cosTheta)); // guard against roundoff past +-1
        arg += std::acos(cosTheta);
        localcenters[bdryList[static_cast<size_t>(k) - 1]] =
            d[static_cast<size_t>(k)] * Complex(std::cos(arg), std::sin(arg));
    }
}

void Packer::updateVdata() {
    const Index layCount = static_cast<Index>(layoutVerts.size());
    if (layCount == 0) {
        std::fprintf(stderr, "Error: \"layoutVerts\" was empty.\n");
    }

    inRadii.assign(static_cast<size_t>(layCount) + 1, {});
    for (Index k = 1; k <= layCount; ++k) {
        Index v = indx2v[k];
        Scalar vrad = localradii[v];
        const auto& flower = flowers[v];
        Index num = vNum[v];
        std::vector<Scalar> data(static_cast<size_t>(num), 0.0);
        Index u = flower[0];
        Scalar urad = localradii[u];
        for (Index j = 1; j <= num; ++j) {
            Scalar wrad = urad;
            u = flower[static_cast<size_t>(j)]; // flower(j+1)
            urad = localradii[u];
            data[static_cast<size_t>(j) - 1] =
                std::sqrt((vrad * urad * wrad) / (vrad + urad + wrad));
        }
        inRadii[k] = std::move(data);
    }

    conduct.assign(static_cast<size_t>(layCount) + 1, 0.0);
    for (Index k = 1; k <= layCount; ++k) {
        Index v = indx2v[k];
        Scalar vrad = localradii[v];
        const auto& iR = inRadii[k];
        Index num = vNum[v];
        const auto& flower = flowers[v];
        Index w = flower[0];
        conduct[k] = (iR[static_cast<size_t>(num) - 1] + iR[0]) / (vrad + localradii[w]);
        for (Index j = 2; j <= num; ++j) {
            Scalar t1 = iR[static_cast<size_t>(j) - 2];
            Scalar t2 = iR[static_cast<size_t>(j) - 1];
            w = flower[static_cast<size_t>(j) - 1];
            conduct[k] += (t1 + t2) / (vrad + localradii[w]);
        }
    }
}

void Packer::layoutCenters() {
    updateVdata();

    const Index layCount = static_cast<Index>(layoutVerts.size());
    const Index colCount = static_cast<Index>(indx2v.size()) - 1 - layCount;

    std::vector<Triplet> tranTriplets;
    tranTriplets.reserve(static_cast<size_t>(tranIJcount));
    for (Index k = 0; k < tranIJcount; ++k) {
        Index row = tranI[static_cast<size_t>(k)];
        Index col = tranJ[static_cast<size_t>(k)];
        Index j = tranJindx[static_cast<size_t>(k)];
        Scalar val;
        if (j > 0) {
            Index v = indx2v[row];
            Scalar vrad = localradii[v];
            Index num = vNum[v];
            const auto& iR = inRadii[v2indx[v]];
            const auto& flower = flowers[v];
            Scalar t1 = (j == 1) ? iR[static_cast<size_t>(num) - 1] : iR[static_cast<size_t>(j) - 2];
            if (j > static_cast<Index>(iR.size())) {
                continue;
            }
            Scalar t2 = iR[static_cast<size_t>(j) - 1];
            Index w = flower[static_cast<size_t>(j) - 1];
            val = ((t1 + t2) / (vrad + localradii[w])) / conduct[v2indx[v]];
        } else {
            val = -1.0;
        }
        tranTriplets.emplace_back(row - 1, col - 1, val);
    }
    transMatrix.resize(layCount, layCount);
    transMatrix.setFromTriplets(tranTriplets.begin(), tranTriplets.end());

    std::vector<Triplet> rhsTriplets;
    rhsTriplets.reserve(static_cast<size_t>(rhsIJcount));
    for (Index k = 0; k < rhsIJcount; ++k) {
        Index row = rhsI[static_cast<size_t>(k)];
        Index col = rhsJ[static_cast<size_t>(k)];
        Index j = rhsJindx[static_cast<size_t>(k)];
        Index v = indx2v[row];
        Scalar vrad = localradii[v];
        Index num = vNum[v];
        const auto& iR = inRadii[v2indx[v]];
        const auto& flower = flowers[v];
        Scalar t1 = (j == 1) ? iR[static_cast<size_t>(num) - 1] : iR[static_cast<size_t>(j) - 2];
        Scalar t2 = iR[static_cast<size_t>(j) - 1];
        Index w = flower[static_cast<size_t>(j) - 1];
        Scalar val = -1.0 * ((t1 + t2) / (vrad + localradii[w])) / conduct[v2indx[v]];
        rhsTriplets.emplace_back(row - 1, col - 1, val);
    }
    rhsMatrix.resize(layCount, colCount);
    rhsMatrix.setFromTriplets(rhsTriplets.begin(), rhsTriplets.end());

    Vector zbReal = Vector::Zero(colCount);
    Vector zbImag = Vector::Zero(colCount);
    for (Index m = 1; m <= colCount; ++m) {
        Complex zb = localcenters[indx2v[layCount + m]];
        zbReal[m - 1] = zb.real();
        zbImag[m - 1] = zb.imag();
    }
    Vector rhsReal = rhsMatrix * zbReal;
    Vector rhsImag = rhsMatrix * zbImag;

    rhs.assign(static_cast<size_t>(layCount) + 1, Complex(0.0, 0.0));
    for (Index k = 1; k <= layCount; ++k) {
        rhs[k] = Complex(rhsReal[k - 1], rhsImag[k - 1]);
    }

    // Solve transMatrix * Z = rhs for real and imaginary parts, reusing one
    // factorization (transMatrix is real, so this is exact and ~2x cheaper
    // than a full complex solve). See SparseLinearSolver.h/.cpp for why
    // Cholesky/CG are never used here: transMatrix is not symmetric.
    SolverStrategy effective = solverStrategy;
    if (effective == SolverStrategy::Auto) {
        effective = (layCount <= 500000) ? SolverStrategy::DirectLU
                                          : SolverStrategy::IterativeBiCGSTAB;
    } else if (effective == SolverStrategy::DirectCholesky ||
               effective == SolverStrategy::IterativeCG) {
        effective = SolverStrategy::DirectLU; // not valid for a non-symmetric system
    }

    Vector zReal(layCount), zImag(layCount);
    bool ok = false;
    if (effective == SolverStrategy::DirectLU) {
        Eigen::SparseLU<SparseMatrix> solver;
        solver.compute(transMatrix);
        if (solver.info() == Eigen::Success) {
            zReal = solver.solve(rhsReal);
            zImag = solver.solve(rhsImag);
            ok = (solver.info() == Eigen::Success);
        }
    } else {
        Eigen::BiCGSTAB<SparseMatrix> solver;
        solver.compute(transMatrix);
        zReal = solver.solve(rhsReal);
        zImag = solver.solve(rhsImag);
        ok = (solver.info() == Eigen::Success);
    }
    if (!ok) {
        throw std::runtime_error(
            "Packer::layoutCenters: sparse solve for interior centers failed to converge");
    }

    for (Index k = 1; k <= layCount; ++k) {
        localcenters[indx2v[k]] = Complex(zReal[k - 1], zImag[k - 1]);
    }
}

void Packer::setEffective() {
    for (size_t k = 0; k < layoutVerts.size(); ++k) {
        Index v = indx2v[static_cast<Index>(k) + 1];
        Scalar targetArea = vAims[v] / 2.0;
        Scalar area = 0.0;
        Index num = vNum[v];
        Complex z = localcenters[v];
        const auto& flower = flowers[v];
        for (Index j = 1; j <= num; ++j) {
            Index jr = flower[static_cast<size_t>(j) - 1];
            Index jl = flower[static_cast<size_t>(j)];
            Complex zr = localcenters[jr];
            Complex zl = localcenters[jl];
            Scalar r = 0.5 * (std::abs(zr - z) + std::abs(zl - z) - std::abs(zr - zl));
            Scalar cC = cosCorner(localcenters[v], localcenters[jr], localcenters[jl]);
            Scalar ang = std::acos(cC);
            area += 0.5 * r * r * ang;
        }
        if (targetArea > 0.001) {
            localradii[v] = std::sqrt(std::max(area / targetArea, 0.0));
        }
    }

    for (Index k = 1; k <= bdryCount; ++k) {
        Index w = bdryList[static_cast<size_t>(k) - 1];
        Scalar targetArea = vAims[w] / 2.0;
        Scalar angsum = 0.0;
        Scalar area = 0.0;
        Index num = vNum[w];
        Complex z = localcenters[w];
        const auto& flower = flowers[w];
        for (Index j = 1; j <= num; ++j) {
            Index jr = flower[static_cast<size_t>(j) - 1];
            Index jl = flower[static_cast<size_t>(j)];
            if (bdryFlags[jr] >= 0 && bdryFlags[jl] >= 0) {
                Complex zr = localcenters[jr];
                Complex zl = localcenters[jl];
                Scalar r = 0.5 * (std::abs(zr - z) + std::abs(zl - z) - std::abs(zr - zl));
                Scalar cC = cosCorner(localcenters[w], localcenters[jr], localcenters[jl]);
                Scalar ang = std::acos(cC);
                angsum += ang;
                area += 0.5 * r * r * ang;
            }
        }

        if (targetArea > 0.001) {
            localradii[w] = std::sqrt(std::max(area / targetArea, 0.0));
        } else if (targetArea < -0.001) {
            localradii[w] = (std::sqrt(std::max(2 * area / angsum, 0.0)) + localradii[w]) / 2.0;
        }
    }
}

Index Packer::continueRiffle(int passNum) {
    constexpr Scalar cutval = 0.01;
    int pass = 0;

    if (passNum <= 0) {
        layoutBdry();
        layoutCenters();
        auto visErr = visualErrors();
        Scalar maxVis = visErr.empty() ? 0.0 : *std::max_element(visErr.begin(), visErr.end());
        visErrMonitor.push_back(maxVis);
        return 0;
    }

    Scalar maxVis = 2 * cutval;
    while (pass < passNum && maxVis > cutval) {
        layoutBdry();
        layoutCenters();
        setEffective();

        auto visErr = visualErrors();
        maxVis = visErr.empty() ? 0.0 : *std::max_element(visErr.begin(), visErr.end());
        visErrMonitor.push_back(maxVis);
        pass++;
    }

    return pass;
}

RiffleResult Packer::riffle(int passNum) {
    RiffleResult result;
    if (mode <= 0) {
        std::fprintf(stderr, "mode = %d suggests there has been an error.\n", mode);
        result.cycles = -1;
        return result;
    }

    auto start = std::chrono::steady_clock::now();
    result.cycles = continueRiffle(passNum);
    auto elapsed = std::chrono::steady_clock::now() - start;
    reapResults();

    if (visErrMonitor.empty()) {
        std::fprintf(stdout, "Riffle: %d passes\n", result.cycles);
    } else {
        result.visMon = visErrMonitor.back();
        result.elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
        std::fprintf(stdout, "\nRiffle: %d passes, visErr %f, elapsed %ld ms\n", result.cycles,
                     result.visMon, result.elapsedMs);
    }
    return result;
}

void Packer::reapResults() {
    centers = localcenters;
    radii = localradii;

    // Spherical packings: recenter so the tangency-point centroid sits at
    // the origin in 3D (the same affine-normalization affineNormalizer/
    // centroid in Geometry.cpp compute), rather than leaving the packing
    // wherever alpha's fixed euclidean-origin placement happens to land
    // once eToSData projects it onto the sphere. This used to happen only
    // inside writepack()'s Spherical branch, which meant it silently
    // depended on going through writepack() to get a well-centered result:
    // any caller reading Packer::centers/radii directly after riffle() --
    // in particular the JNI bridge's computeMaximalPackingFromComplex,
    // which never calls writepack() at all -- got the un-normalized,
    // possibly lopsided placement instead. Doing it here means every
    // caller gets a normalized packing for free. writepack() (PackerIO.cpp)
    // deliberately still does its own affineNormalizer call too, as a
    // defensive no-op for a caller that writes a Spherical packing without
    // riffling it first (so this step never ran) -- affineNormalizer is
    // idempotent on an already-centered input (its very first check is
    // already within tolerance, so it returns the identity transform
    // immediately), so that redundancy costs one cheap extra pass on the
    // normal riffle-then-write path, not real duplicated work. Hyperbolic
    // has no equivalent step: writepack()'s Hyperbolic branch is a plain
    // per-vertex eToHData conversion with no affine pre-step to move.
    if (hes == Geometry::Spherical) {
        std::vector<Complex> T = loadTangency(centers, radii);
        auto [A, B] = geom::affineNormalizer(T);
        for (Index v = 1; v <= nodeCount; ++v) {
            centers[v] = A * centers[v] + B;
            radii[v] = A * radii[v];
        }
    }
}

std::vector<Scalar> Packer::visualErrors() const {
    std::vector<Scalar> visualErr(layoutVerts.size(), 0.0);
    for (size_t k = 0; k < layoutVerts.size(); ++k) {
        Index v = layoutVerts[k];
        Complex centv = localcenters[v];
        Scalar radv = localradii[v];
        const auto& flower = flowers[v];
        Index num = vNum[v];
        Scalar maxerr = 0.0;
        for (Index j = 1; j <= num; ++j) {
            Index w = flower[static_cast<size_t>(j) - 1];
            Scalar cdiff = std::abs(centv - localcenters[w]);
            Scalar raddiff = radv + localradii[w];
            Scalar me = std::abs(cdiff - raddiff) / radv;
            if (me > maxerr) maxerr = me;
        }
        visualErr[k] = maxerr;
    }
    return visualErr;
}

std::vector<Scalar> Packer::angsumErrorsImpl(const std::vector<Scalar>* radiiOverride) const {
    const std::vector<Scalar>& r = radiiOverride ? *radiiOverride : localradii;
    std::vector<Scalar> diffs(layoutVerts.size(), 0.0);
    for (size_t k = 0; k < layoutVerts.size(); ++k) {
        Index v = indx2v[static_cast<Index>(k) + 1];
        Scalar diff = -2.0 * kPi;
        Scalar rv = r[v];
        Index num = vNum[v];
        const auto& flower = flowers[v];
        Index u = flower[0];
        for (Index j = 1; j <= num; ++j) {
            Index w = u;
            u = flower[static_cast<size_t>(j)];
            diff += std::acos(cosAngle(rv, r[w], r[u]));
        }
        diffs[k] = diff;
    }
    return diffs;
}

void Packer::packStatus() const {
    int hesIdx = static_cast<int>(hes) + 1; // -1,0,1 -> 0,1,2
    std::fprintf(stdout, "Status for %s packing \"%s\":\n", kHesNames[hesIdx], fileName.c_str());
    std::fprintf(stdout, "  Mode is \"%s\"", kPackModes[mode - 1]);
    if (mode == 2) {
        std::fprintf(stdout, ", with corner vertices ");
        for (Index c : corners) std::fprintf(stdout, "%d ", c);
    }
    std::fprintf(stdout, ".\n");
    std::fprintf(stdout,
                 "  Packing has %d vertices, of which %zu are subject to packing adjustments.\n",
                 nodeCount, layoutVerts.size());
    auto verr = visualErrors();
    Scalar mx = verr.empty() ? 0.0 : *std::max_element(verr.begin(), verr.end());
    std::fprintf(stdout, "  Maximum visual error is %f. For details see \"visualErrors\" array.\n",
                 mx);
}

} // namespace gopack
