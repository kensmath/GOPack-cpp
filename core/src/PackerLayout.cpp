#include "gopack/Packer.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
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
    if (m < 1 || m > 2) {
        std::fprintf(stderr, "Mode choices: 1 = max packing, 2 = polygonal packing.\n");
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

    (void)crns;
    (void)angs;
    mode = -1;
    throw NotImplementedError(
        "Packer::setMode: polygonal packing mode (mode 2) has not been ported "
        "yet (see setMode.m's mode==2 branch, and setPolyCenters.m / "
        "setRectCenters.m). The scoped first pass of this port covers "
        "maximal packing (mode 1) only.");
}

void Packer::layoutBdry() {
    if (mode == 2) {
        setPolyCenters();
        return;
    }
    setHoroCenters();
}

void Packer::setRectCenters() {
    throw NotImplementedError(
        "Packer::setRectCenters: not yet ported (setRectCenters.m); "
        "polygonal/rectangle packing mode is deferred.");
}

void Packer::setPolyCenters() {
    throw NotImplementedError(
        "Packer::setPolyCenters: not yet ported (setPolyCenters.m); "
        "polygonal/rectangle packing mode is deferred.");
}

Scalar Packer::getAspect() {
    throw NotImplementedError(
        "Packer::getAspect: not yet ported (getAspect.m); only meaningful for "
        "the not-yet-ported polygonal/rectangle packing mode.");
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
