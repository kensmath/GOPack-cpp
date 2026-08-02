#pragma once
//
// GOPack C++ -- thin wrapper around the sparse linear solve backend.
//
// This isolates the rest of the codebase from the specific sparse library in
// use. Today that's Eigen's built-in sparse solvers (SimplicialLDLT for
// symmetric positive-definite systems, SparseLU otherwise, BiCGSTAB as an
// iterative fallback for very large systems where a direct factorization is
// too memory-hungry). If GOPACK_USE_SUITESPARSE is enabled at configure time,
// this is where the CHOLMOD/UMFPACK backend gets swapped in instead -- the
// public API below does not change.
//
#include "gopack/Types.h"

namespace gopack {

enum class SolverStrategy {
    Auto,          // pick based on matrix size/structure
    DirectCholesky,
    DirectLU,
    IterativeCG,
    IterativeBiCGSTAB
};

class SparseLinearSolver {
public:
    explicit SparseLinearSolver(SolverStrategy strategy = SolverStrategy::Auto);

    // Solves A x = b. Returns false (and leaves x unchanged) if the solver
    // fails to converge / factorize.
    bool solve(const SparseMatrix& A, const Vector& b, Vector& x);

    // Diagnostics from the most recent solve() call.
    int iterations() const { return lastIterations_; }
    Scalar residualNorm() const { return lastResidualNorm_; }

private:
    SolverStrategy strategy_;
    int lastIterations_ = 0;
    Scalar lastResidualNorm_ = 0.0;
};

} // namespace gopack
