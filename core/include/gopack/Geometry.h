#pragma once
//
// GOPack C++ -- geometry conversion utilities.
//
// Direct translations of the free-standing GOPack MATLAB helper functions
// (code/cosAngle.m, cosCorner.m, e_to_h_data.m, h_to_e_data.m, e_to_s_data.m,
// s_to_e_data.m, sph_tangent.m, proj_vec_to_s.m, s_pt_to_vec.m, Centroid.m,
// affineNormalizer.m). GOPack always computes in euclidean geometry
// internally; these are used only when reading/writing hyperbolic or
// spherical packings.
#include <array>
#include <vector>

#include "gopack/Types.h"

namespace gopack::geom {

// cosAngle.m -- cosine of the angle at a circle of radius r in a triple of
// mutually tangent circles with radii r, r1, r2.
Scalar cosAngle(Scalar r, Scalar r1, Scalar r2);

// cosCorner.m -- cosine of the angle at z1 in the triangle (z1, z2, z3).
Scalar cosCorner(Complex z1, Complex z2, Complex z3);

// e_to_h_data.m / h_to_e_data.m -- euclidean <-> hyperbolic (Poincare disc)
// circle data. Hyperbolic radius is negative for horocycles, by convention,
// with -hr equal to the euclidean radius in that case.
std::pair<Complex, Scalar> eToHData(Complex ez, Scalar er);
std::pair<Complex, Scalar> hToEData(Complex hz, Scalar hr);

// e_to_s_data.m / s_to_e_data.m -- euclidean <-> spherical circle data.
// Spherical centers are encoded as complex (theta, phi) polar coordinates.
std::pair<Complex, Scalar> eToSData(Complex ez, Scalar er);
std::pair<Complex, Scalar> sToEData(Complex sz, Scalar sr);

// sph_tangent.m -- unit tangent vector at ctr1 pointing toward ctr2 on S^2.
std::array<Scalar, 3> sphTangent(Complex ctr1, Complex ctr2);

// proj_vec_to_s.m / s_pt_to_vec.m -- project a 3-vector onto S^2 (returned
// in (theta,phi) polar form) and the inverse conversion.
Complex projVecToS(const std::array<Scalar, 3>& vec);
std::array<Scalar, 3> sPtToVec(Complex sz);

struct CentroidResult {
    Scalar normSq = 0.0;
    Scalar x = 0.0, y = 0.0, z = 0.0;
};

// Centroid.m -- centroid (in 3-space, via stereographic projection) of
// points P after applying the affine map z -> trans[0]*z + trans[1] + i*trans[2].
CentroidResult centroid(const std::vector<Complex>& P, const std::array<Scalar, 3>& trans);

// affineNormalizer.m -- find real a, complex b such that the affine map
// z -> a*z + b moves the centroid of T (after stereographic projection) to
// the origin. Used to normalize spherical output. T is passed by value
// because the search mutates its own working copy.
std::pair<Scalar, Complex> affineNormalizer(std::vector<Complex> T);

} // namespace gopack::geom
