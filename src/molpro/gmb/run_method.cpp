#include "run_method.h"
#include "problem_eom-ccsd.h"

#include "problem_ccsd.h"
#include "problem_ccd.h"
#include "problem_mp2.h"

#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory.h>

#include <molpro/iostream.h>

void run_gs(hamiltonian<> &ham, const std::string &method, std::unique_ptr<problem_gen> &problem, std::unique_ptr<amplitudes<>> &ptampl) {

  // set CCSD amplitudes
  if (method.find("ccs") != std::string::npos)
    ptampl->set(t1, container(ham.m2get(f_ov).get_space()));
  if ((method.find("ccsd") != std::string::npos) || (method.find("ccd") != std::string::npos))
    ptampl->set(t2, container(ham.m4get(i_oovv).get_space()));
  else if (method.find("mp2") != std::string::npos)
    ptampl->set(t2, ham.m4get(i_oovv));

  // set problem
  if (method.find("ccsd") != std::string::npos)
    problem = std::make_unique<problem_ccsd>(ham);
  else if (method.find("ccd") != std::string::npos)
    problem = std::make_unique<problem_ccd>(ham);
  else if (method.find("mp2") != std::string::npos)
    problem = std::make_unique<problem_mp2>(ham);

  molpro::cout << "Running " << *problem << std::endl;

  // set solver
  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<amplitudes<>>("DIIS", "max_size_qspace=8");
  auto residual = *ptampl;

  // solver options
  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Iteration);
  // solver->set_convergence_threshold(1.0e-14);

  solver->solve(*ptampl, residual, *problem);
  solver->solution(*ptampl, residual);
  problem->energy(*ptampl);
}

void run_es(const hamiltonian<> &ham, const std::string &method, std::unique_ptr<problem_eom> &problem, const std::unique_ptr<amplitudes<>> &ptampl, const size_t &nroots) {
  molpro::cout << "\nRunning EOM-CCSD" << std::endl;
  
  // set EOM-CCSD amplitudes
  problem = std::make_unique<problem_eom_ccsd>(ham, *ptampl);

  std::unique_ptr<amplitudes<>> prampl{std::make_unique<amplitudes<>>()};
  prampl->set(r1, container(ptampl->m2get(t1).get_space()));
  prampl->set(r2, container(ptampl->m4get(t2).get_space()));
  std::vector<amplitudes<>> v_rampl(nroots, *prampl);
  
  // set solver 
  auto solver = molpro::linalg::itsolv::create_LinearEigensystem<amplitudes<>>("Davidson");
  auto residuals_es = v_rampl;
  
  // set options
  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Iteration);
  solver->set_n_roots(nroots);
  solver->set_convergence_threshold(1.0e-5);

  // solve
  solver->solve(v_rampl, residuals_es, *problem, true);
  problem->set_energy(solver->eigenvalues());
  std::vector<int> v_nroots(nroots);
  for (size_t i = 0; i < nroots; i++)
    v_nroots[i] = i;
  solver->solution(v_nroots, v_rampl, residuals_es);
  problem->character(v_rampl);

  // check for zero elements
  auto energies = problem->get_energy();

  for (size_t i = 0; i < energies.size(); i++) {
    if (energies[i] < 10e-5) 
    {
      std::cout << "\nWarning: Found 0 eigenvalue!" << std::endl;
      problem->check_eigenvalue(v_rampl[i]);
      
      // set EOM-CCSD amplitudes
      // problem = std::make_unique<problem_eom_ccsd>(ham, *ptampl);

      // std::unique_ptr<amplitudes<>> prampl{std::make_unique<amplitudes<>>()};
      prampl->set(r1, container(v_rampl[i].m2get(r1)));
      prampl->set(r2, container(v_rampl[i].m4get(r2)));
      std::vector<amplitudes<>> v_rampl(1, *prampl);

      // set solver 
      auto solver = molpro::linalg::itsolv::create_LinearEigensystem<amplitudes<>>("Davidson");
      auto residuals_es = v_rampl;

      // set options
      solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Detailed);
      solver->set_n_roots(1);
      // solver->set_max_iter(1);
      solver->set_convergence_threshold(1.0e-10);

      // solve
      solver->solve(v_rampl, residuals_es, *problem, false);
      problem->set_energy(solver->eigenvalues());

    }
  }
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<amplitudes<>, amplitudes<>>;