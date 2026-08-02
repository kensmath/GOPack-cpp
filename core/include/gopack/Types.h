#pragma once
//
// GOPack C++ -- core scalar/index/geometry type definitions.
//
#include <complex>
#include <cstdint>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace gopack {

using Index   = int;
using Scalar  = double;
using Complex = std::complex<double>;

using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::ColMajor, Index>;
using Triplet      = Eigen::Triplet<Scalar, Index>;

// Mirrors GOPacker's `hes` field: -1 hyperbolic, 0 euclidean, +1 spherical.
enum class Geometry : int {
    Hyperbolic = -1,
    Euclidean  = 0,
    Spherical  = 1
};

constexpr double kPi    = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;

} // namespace gopack
