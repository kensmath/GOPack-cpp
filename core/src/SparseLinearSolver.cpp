#include "gopack/SparseLinearSolver.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

namespace gopack {

SparseLinearSolver::SparseLinearSolver(SolverStrategy strategy) : strategy_(strategy) {}

bool SparseLinearSolver::solve(const SparseMatrix& A, const Vector& b, Vector& x) {
    lastIterations_ = 0;
    lastResidualNorm_ = 0.0;

    auto strategy = strategy_;
    if (strategy == SolverStrategy::Auto) {
        // IMPORTANT: the system GOPack's layoutCenters() assembles
        // (transMatrix) is a *non-symmetric* diagonally-dominant operator --
        // diagonal -1, row entries for interior petals only (boundary petals
        // are moved to the right-hand side), so it is not, in general,
        // symmetric positive definite. ConjugateGradient/Cholesky both
        // silently assume symmetry (Eigen's SimplicialLDLT only reads one
        // triangle of A), so defaulting to either here would be a correctness
        // bug, not just a performance choice. Auto therefore only ever
        // chooses between a general direct LU factorization (matches
        // MATLAB's mldivide behavior on this matrix most closely) and a
        // general iterative method (BiCGSTAB) for very large systems where a
        // direct factorization would be too memory-hungry.
        strategy = (A.rows() <= 500000) ? SolverStrategy::DirectLU
                                         : SolverStrategy::IterativeBiCGSTAB;
    }

    switch (strategy) {
        case SolverStrategy::DirectCholesky: {
            Eigen::SimplicialLDLT<SparseMatrix> solver;
            solver.compute(A);
            if (solver.info() != Eigen::Success) {
                return false;
            }
            x = solver.solve(b);
            lastResidualNorm_ = (A * x - b).norm();
            return solver.info() == Eigen::Success;
        }
        case SolverStrategy::DirectLU: {
            Eigen::SparseLU<SparseMatrix> solver;
            solver.compute(A);
            if (solver.info() != Eigen::Success) {
                return false;
            }
            x = solver.solve(b);
            lastResidualNorm_ = (A * x - b).norm();
            return solver.info() == Eigen::Success;
        }
        case SolverStrategy::IterativeCG: {
            Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> solver;
            solver.compute(A);
            x = solver.solve(b);
            lastIterations_ = static_cast<int>(solver.iterations());
            lastResidualNorm_ = solver.error();
            return solver.info() == Eigen::Success;
        }
        case SolverStrategy::IterativeBiCGSTAB: {
            Eigen::BiCGSTAB<SparseMatrix> solver;
            solver.compute(A);
            x = solver.solve(b);
            lastIterations_ = static_cast<int>(solver.iterations());
            lastResidualNorm_ = solver.error();
            return solver.info() == Eigen::Success;
        }
        default:
            return false;
    }
}

} // namespace gopack
