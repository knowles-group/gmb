#include "problem_ccsd.h"
#include <iostream>
#include <molpro/linalg/itsolv/SolverFactory.h>

int main(int argc, char* argv[]) {
  supercontainer ham;
  auto problem = problem_ccsd(ham);
  using Rvector = problem_ccsd::container_t;
  using Qvector = Rvector;
  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<Rvector, Qvector>("DIIS");
  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Summary);
  Rvector c, g;
  if (not solver->solve(c, g, problem))
    std::cout << "failed\n";
  else
    std::cout << "converged in " << solver->statistics().iterations << " iterations\n";
  solver->solution(c, g);
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<supercontainer<>, supercontainer<>, std::map<size_t, double>>;
