#include "gopack/Geometry.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

namespace gopack::geom {

Scalar cosAngle(Scalar r, Scalar r1, Scalar r2) {
    Scalar c = r1 * r2;
    return 1.0 - 2.0 * c / (r * r + r * (r1 + r2) + c);
}

Scalar cosCorner(Complex z1, Complex z2, Complex z3) {
    Scalar l2 = std::abs(z2 - z1);
    Scalar l3 = std::abs(z3 - z1);
    Scalar l23 = std::abs(z3 - z2);
    Scalar denom = 2.0 * l2 * l3;
    Scalar cosang = (l2 * l2 + l3 * l3 - l23 * l23) / denom;
    return std::max(std::min(cosang, 1.0), -1.0);
}

std::pair<Complex, Scalar> eToHData(Complex ez, Scalar er) {
    Scalar aec = std::abs(ez);
    Scalar dist = aec + er;
    if (dist > 1.000000000001) { // not in closed disc; push in to horocycle
        aec = aec / dist;
        er = er / dist;
        dist = 1.0;
    }

    // is this a horocycle?
    if (0.99999999 < dist) {
        if (std::abs(er) > 0.99999999) {
            er = er / 2.0;
        }
        Scalar hr = -er;
        Complex hz;
        if (aec < 0.0001) {
            hz = Complex(0.0, 0.0);
        } else {
            hz = ez * (1.0 / aec);
        }
        return {hz, hr};
    }

    Scalar c2 = aec * aec;
    Scalar r2 = er * er;
    Complex hz;
    if (aec < 0.0000000000001) {
        hz = Complex(0.0, 0.0);
    } else {
        Scalar t = 1 + c2 - r2;
        Scalar b = std::sqrt((t + 2 * aec) / (t - 2 * aec));
        Scalar ahc = (b - 1) / (b + 1);
        hz = ez * (ahc / aec);
    }

    Scalar t = 1 + r2 - c2;
    Scalar s = std::sqrt((t - 2 * er) / (t + 2 * er)); // s-radius
    Scalar x = 1.0 - s * s;                            // x-radius
    Scalar hr;
    if (x > 0.0001) {
        hr = (-0.5) * std::log(1.0 - x);
    } else {
        hr = x * (1.0 + x * (0.5 + x / 3)) / 2; // 2nd order approximation
    }
    return {hz, hr};
}

std::pair<Complex, Scalar> hToEData(Complex hz, Scalar hr) {
    if (hr < 0) { // horocycle
        Scalar er = -hr;
        Complex ez = hz * (1.0 - er);
        return {ez, er};
    }

    Scalar ahc = std::abs(hz);
    Scalar sRad = std::exp(-hr);
    Scalar n1 = (1 + sRad) * (1 + sRad);
    Scalar n2 = n1 - ahc * ahc * (1 - sRad) * (1 - sRad);
    Scalar er = (1.0 - sRad * sRad) * (1.0 - ahc * ahc) / n2;
    Scalar b = 4.0 * sRad / n2;
    Complex ez = hz * b;
    return {ez, er};
}

std::array<Scalar, 3> sPtToVec(Complex sz) {
    Scalar s = std::sin(sz.imag());
    Scalar c = std::cos(sz.real());
    return {s * c, s * std::sin(sz.real()), std::cos(sz.imag())};
}

Complex projVecToS(const std::array<Scalar, 3>& vec) {
    constexpr Scalar kTol = 0.0000000000001;
    Scalar dist = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    if (dist < kTol) {
        return Complex(0.0, 0.0);
    }
    return Complex(std::atan2(vec[1], vec[0]), std::acos(vec[2] / dist));
}

std::array<Scalar, 3> sphTangent(Complex ctr1, Complex ctr2) {
    constexpr Scalar kTol = 0.00000000001;
    std::array<Scalar, 3> A = sPtToVec(ctr1);
    std::array<Scalar, 3> B = sPtToVec(ctr2);
    Scalar d = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
    std::array<Scalar, 3> P = {B[0] - d * A[0], B[1] - d * A[1], B[2] - d * A[2]};

    Scalar vn = std::sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
    if (vn < kTol) {
        Scalar pn = std::sqrt(A[1] * A[1] + A[2] * A[2]);
        if (pn > 0.001) {
            return {0.0, A[1] / pn, -1.0 * A[2] / pn};
        }
        return {1.0, 0.0, 0.0};
    }
    return {P[0] / vn, P[1] / vn, P[2] / vn};
}

std::pair<Complex, Scalar> eToSData(Complex ez, Scalar er) {
    constexpr Scalar kTol = 0.00000000001;
    Scalar ns = ez.real() * ez.real() + ez.imag() * ez.imag();
    Scalar rr = std::abs(er);

    if (rr < kTol) {
        // er too small; project center, er unchanged.
        Scalar denom = ns + 1.0;
        Scalar tmpd = 1.0 / denom;
        std::array<Scalar, 3> P3 = {(2 * ez.real()) * tmpd, (2 * ez.imag()) * tmpd,
                                     (2.0 - denom) * tmpd};
        if (P3[2] > (1.0 - kTol)) {
            return {Complex(0.0, 0.0), er};
        }
        if (P3[2] < (kTol - 1.0)) {
            return {Complex(0.0, kPi), er};
        }
        return {Complex(std::atan2(P3[1], P3[0]), std::acos(P3[2])), er};
    }

    Scalar norm = std::sqrt(ns);
    Scalar x, y, a, b, mn;
    if (norm < kTol) {
        mn = -rr;
        x = mn;
        y = 0.0;
        a = rr;
        b = 0.0;
    } else {
        Scalar denom = 1 / norm;
        mn = norm - rr;
        x = mn * ez.real() * denom;
        y = mn * ez.imag() * denom;
        a = (norm + rr) * ez.real() * denom;
        b = (norm + rr) * ez.imag() * denom;
    }

    Scalar d1 = (x * x + y * y + 1.0);
    Scalar tmpd = 1.0 / d1;
    std::array<Scalar, 3> P1 = {2.0 * x * tmpd, 2.0 * y * tmpd, (2.0 - d1) * tmpd};
    Scalar d2 = a * a + b * b + 1.0;
    tmpd = 1.0 / d2;
    std::array<Scalar, 3> P2 = {2.0 * a * tmpd, 2.0 * b * tmpd, (2.0 - d2) * tmpd};

    constexpr Scalar kBrk = 100.0 * kTol;
    bool midflag = false;
    std::array<Scalar, 3> P3 = {0, 0, 0};
    if (mn <= -kBrk) {
        midflag = true;
        P3 = {0.0, 0.0, 1.0};
    } else if (mn <= kBrk && norm > 2) {
        midflag = true;
        P3 = {ez.real() / norm, ez.imag() / norm, 0.0};
    }

    Scalar rad;
    std::array<Scalar, 3> E{};
    if (midflag) {
        Scalar dd1 = P1[0] * P3[0] + P1[1] * P3[1] + P1[2] * P3[2];
        if (dd1 >= 1.0) dd1 = 1.0 - kTol;
        Scalar dd2 = P2[0] * P3[0] + P2[1] * P3[1] + P2[2] * P3[2];
        if (dd2 >= 1.0) dd2 = 1.0 - kTol;
        Scalar ang13 = std::acos(dd1);
        Scalar ang23 = std::acos(dd2);
        rad = (ang13 + ang23) / 2.0;
        E = (ang13 < ang23) ? P1 : P2;
        Complex v(std::atan2(E[1], E[0]), std::acos(E[2]));
        Complex w(std::atan2(P3[1], P3[0]), std::acos(P3[2]));
        // T computed below (shared code path)
        auto T = sphTangent(v, w);
        std::array<Scalar, 3> C = {E[0] * std::cos(rad) + T[0] * std::sin(rad),
                                    E[1] * std::cos(rad) + T[1] * std::sin(rad),
                                    E[2] * std::cos(rad) + T[2] * std::sin(rad)};
        Scalar sr = rad;
        if (rad < 0) {
            sr = kPi - rad;
            C = {-1.0, -1.0, -1.0};
        }
        Complex sz;
        if (C[2] > 1 - kTol) {
            sz = Complex(0.0, 0.0);
        } else if (C[2] < (kTol - 1.0)) {
            sz = Complex(0.0, kPi);
        } else {
            sz = Complex(std::atan2(C[1], C[0]), std::acos(C[2]));
        }
        return {sz, sr};
    }

    Scalar dd1 = P1[0] * P2[0] + P1[1] * P2[1] + P1[2] * P2[2];
    if (dd1 >= 1.0) dd1 = 1.0 - kTol;
    rad = std::acos(dd1) / 2.0;
    E = P1;
    Complex v(std::atan2(E[1], E[0]), std::acos(E[2]));
    Complex w(std::atan2(P2[1], P2[0]), std::acos(P2[2]));
    auto T = sphTangent(v, w);
    std::array<Scalar, 3> C = {E[0] * std::cos(rad) + T[0] * std::sin(rad),
                                E[1] * std::cos(rad) + T[1] * std::sin(rad),
                                E[2] * std::cos(rad) + T[2] * std::sin(rad)};
    Scalar sr = rad;
    if (rad < 0) {
        sr = kPi - rad;
        C = {-1.0, -1.0, -1.0};
    }
    Complex sz;
    if (C[2] > 1 - kTol) {
        sz = Complex(0.0, 0.0);
    } else if (C[2] < (kTol - 1.0)) {
        sz = Complex(0.0, kPi);
    } else {
        sz = Complex(std::atan2(C[1], C[0]), std::acos(C[2]));
    }
    return {sz, sr};
}

std::pair<Complex, Scalar> sToEData(Complex sz, Scalar sr) {
    constexpr Scalar kTol = 0.0000000000001;
    Scalar flipflag = 1.0;
    std::array<Scalar, 3> V = sPtToVec(sz);
    Scalar phi = sz.imag();

    // essentially hits infinity (south pole)?
    if (std::abs(phi + sr - kPi) < kTol) {
        sr = sr + 2.0 * kTol;
    }

    // encloses infinity?
    if ((phi + sr) >= (kPi + kTol)) {
        sr = kPi - sr;
        V = {-V[0], -V[1], -V[2]};
        sz = geom::projVecToS(V);
        flipflag = -1.0;
    }

    Scalar up = phi + sr;
    Scalar down = phi - sr;

    // essentially centered at north pole?
    if (std::abs(phi) < kTol) {
        Scalar er = std::sin(up) / (1.0 + std::cos(up));
        Complex ez(0.0, 0.0);
        if (flipflag < 0) er = er * -1.0;
        return {ez, er};
    }

    // essentially centered at south pole
    if (std::abs(phi - kPi) < kTol) {
        return {Complex(0.0, 0.0), -100000.0};
    }

    // circle essentially passes through infinity; decrease 'up' slightly
    if (std::abs(up - kPi) < 0.00001) {
        up = up - 0.00001;
    }

    Scalar RR = std::sin(up) / (1.0 + std::cos(up));
    Scalar rr = std::sin(down) / (1.0 + std::cos(down));
    Scalar er = std::abs(RR - rr) / 2.0;
    if (flipflag < 0) {
        er = er * -1.0;
        Scalar m = (RR + rr) / 2.0;
        Complex ez(V[0] * m / std::sin(phi), V[1] * m / std::sin(phi));
        return {ez, er};
    }
    return {Complex(0.0, 0.0), er};
}

CentroidResult centroid(const std::vector<Complex>& P, const std::array<Scalar, 3>& trans) {
    CentroidResult res;
    Scalar sumX = 0, sumY = 0, sumZ = 0;
    const size_t n = P.size();
    for (const auto& p : P) {
        Scalar mu = trans[0] * p.real() + trans[1];
        Scalar mv = trans[0] * p.imag() + trans[2];
        Scalar sq = mu * mu + mv * mv;
        Scalar denom = 1 + sq;
        sumX += 2 * mu / denom;
        sumY += 2 * mv / denom;
        sumZ += (1 - sq) / denom;
    }
    res.x = sumX / static_cast<Scalar>(n);
    res.y = sumY / static_cast<Scalar>(n);
    res.z = sumZ / static_cast<Scalar>(n);
    res.normSq = res.x * res.x + res.y * res.y + res.z * res.z;
    return res;
}

std::pair<Scalar, Complex> affineNormalizer(std::vector<Complex> T) {
    std::array<Scalar, 3> M = {1.0, 0.0, 0.0};
    Scalar bestsq = centroid(T, M).normSq;
    constexpr Scalar kNTol = 0.001;
    constexpr int kCycles = 20;

    // Deviation from affineNormalizer.m: floor on the scale coordinate
    // (m[0]/j==0 below), strictly greater than zero. centroid()'s objective
    // only ever uses trans[0] *squared* (mu=trans[0]*p.real()+trans[1], sq=
    // mu*mu+mv*mv), so it cannot distinguish a scale of 'a' from '-a' --
    // nothing in the literal source stops this greedy coordinate search from
    // wandering across that sign boundary, and it's not just a hypothetical
    // tie: direct stress-testing (200k random point clouds, see
    // HANDOFFrandomtrinorepack.md's follow-up investigation) shows roughly
    // 0.4% of inputs make the unguarded search land on a genuinely negative
    // scale. That matters here specifically because the sole caller,
    // Packer::reapResults(), applies this scale uniformly to every vertex's
    // radius ("radii[v] = A * radii[v]") -- a negative (or zero) scale
    // silently flips every radius negative (or collapses them all to zero)
    // for the whole packing, which is exactly the failure
    // tests/test_sphere_normalize.cpp caught in CI. kMinScale is set far
    // below any legitimate correction this near-identity recentering step
    // should ever need (it only nudges an already-good packing), so this
    // only ever blocks the pathological sign-flip/collapse case, never a
    // real in-range correction.
    constexpr Scalar kMinScale = 1e-6;

    int outercount = 0;
    while (bestsq > kNTol && outercount < kCycles) {
        Scalar delt = 2.0;
        std::array<Scalar, 3> m = {1.0, 0.0, 0.0};
        int count = 0;

        while (bestsq > kNTol && count < kCycles) {
            int gotOne = 0;
            for (int j = 0; j < 3; ++j) {
                Scalar holdp = m[j];
                m[j] = m[j] + delt;
                Scalar newnorm = centroid(T, m).normSq;
                m[j] = holdp;
                if (newnorm < bestsq) {
                    bestsq = newnorm;
                    gotOne = j + 1;
                } else if (j != 0 || holdp - delt > kMinScale) {
                    m[j] = m[j] - delt;
                    newnorm = centroid(T, m).normSq;
                    m[j] = holdp;
                    if (newnorm < bestsq) {
                        bestsq = newnorm;
                        gotOne = -(j + 1);
                    }
                }
            }

            if (gotOne == 0) {
                delt = delt / 2;
            } else {
                int idx = std::abs(gotOne) - 1;
                m[idx] += (gotOne > 0) ? delt : -delt;
            }
            ++count;
        }

        if (bestsq < kNTol) {
            std::array<Scalar, 3> newM = {m[0] * M[0], m[0] * M[1] + m[1], m[0] * M[2] + m[2]};
            M = newM;
            return {M[0], Complex(M[1], M[2])};
        } else {
            for (auto& v : T) {
                v = m[0] * v + Complex(m[1], m[2]);
            }
            std::array<Scalar, 3> newM = {m[0] * M[0], m[0] * M[1] + m[1], m[0] * M[2] + m[2]};
            M = newM;
        }
        ++outercount;
    }

    // Defensive backstop, not expected to ever trigger given the per-step
    // guard above (each outer round's m[0] is kept > kMinScale throughout,
    // so the accumulated product M[0] -- a product of positive factors --
    // must itself stay positive): if M[0] somehow isn't strictly positive
    // anyway, fall back to the identity transform rather than ever handing
    // reapResults() a scale that would corrupt every radius.
    if (!(M[0] > 0.0)) {
        std::fprintf(stderr,
                     "affineNormalizer: internal scale guard tripped (M[0]=%g) -- "
                     "returning identity transform instead\n",
                     M[0]);
        return {1.0, Complex(0.0, 0.0)};
    }
    return {M[0], Complex(M[1], M[2])};
}

std::vector<Complex> randBdryPts(std::vector<Complex> graph, Index M) {
    Index graphCount = static_cast<Index>(graph.size());
    if (graphCount < 3 || M < 3) {
        std::fprintf(stderr, "Poor data in \"randBdryPts\"\n");
        return {};
    }

    // close up if necessary
    if (std::abs(graph.front().real() - graph.back().real()) > 0.001 &&
        std::abs(graph.front().imag() - graph.back().imag()) > 0.001) {
        graph.push_back(graph.front());
    }
    graphCount = static_cast<Index>(graph.size());

    // mark off by polygon (arc) length
    std::vector<Scalar> lengthMarks(static_cast<size_t>(graphCount), 0.0);
    for (Index i = 1; i < graphCount; ++i) {
        lengthMarks[static_cast<size_t>(i)] =
            lengthMarks[static_cast<size_t>(i) - 1] +
            std::abs(graph[static_cast<size_t>(i)] - graph[static_cast<size_t>(i) - 1]);
    }

    // find M random ordered param spots in [0, total length]
    static thread_local std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<Scalar> unit(0.0, 1.0);
    std::vector<Scalar> arcSpots(static_cast<size_t>(M));
    for (Index i = 0; i < M; ++i) {
        arcSpots[static_cast<size_t>(i)] = unit(rng) * lengthMarks.back();
    }
    std::sort(arcSpots.begin(), arcSpots.end());

    // convert by interpolation to points on the graph
    std::vector<Complex> result(static_cast<size_t>(M));
    Index spot = 0; // segment [spot, spot+1], 0-indexed into graph/lengthMarks
    Scalar lastLength = lengthMarks[static_cast<size_t>(spot)];
    Scalar nextLength = lengthMarks[static_cast<size_t>(spot) + 1];
    for (Index i = 0; i < M; ++i) {
        while (arcSpots[static_cast<size_t>(i)] < lastLength && spot > 0) {
            spot--;
            lastLength = lengthMarks[static_cast<size_t>(spot)];
            nextLength = lengthMarks[static_cast<size_t>(spot) + 1];
        }
        // Defensive addition beyond the literal source (see header comment):
        // guard against ever indexing lengthMarks[spot+1] out of bounds in
        // the (essentially unreachable in practice, since arcSpots is
        // strictly < lengthMarks.back() almost surely) edge case where an
        // arc spot lands exactly on the total path length.
        while (arcSpots[static_cast<size_t>(i)] > nextLength && spot + 2 < graphCount) {
            spot++;
            lastLength = lengthMarks[static_cast<size_t>(spot)];
            nextLength = lengthMarks[static_cast<size_t>(spot) + 1];
        }

        Scalar setlength = nextLength - lastLength;
        Scalar ratio = (setlength > 0.0) ? (arcSpots[static_cast<size_t>(i)] - lastLength) / setlength : 0.0;
        result[static_cast<size_t>(i)] =
            graph[static_cast<size_t>(spot)] +
            ratio * (graph[static_cast<size_t>(spot) + 1] - graph[static_cast<size_t>(spot)]);
    }
    return result;
}

} // namespace gopack::geom
