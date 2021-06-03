#include "init.h"
#include "problem_ccsd.h"
#include "problem_ccd.h"
#include "problem_eom-ccsd.h"
#include "hamiltonian.h"
#include "amplitudes.h"
#include "expressions/energy_hf.h"
#include <memory>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <chrono>
#include <ctime>

std::string filename;
// test case
// std::string test_case = "hubbard"; // electron only
std::string test_case = "polariton"; // electron + photons
// std::string test_case = "smallcoupling_newd";
// std::string test_case = "zerocoupling";

using namespace gmb;
int main(int argc, char const *argv[]) {
  std::ios_base::sync_with_stdio(false);
  filename = argv[0];
  if (filename.find_last_of("/") != std::string::npos)
    filename.resize(filename.find_last_of("/"));
  filename += "/"+test_case+".fcidump";
  auto start = std::chrono::system_clock::now();

  hamiltonian<> hamiltonian;
  std::unique_ptr<amplitudes<>> ptampl(new amplitudes());
  std::unique_ptr<problem_gen> problem;
  std::string method_gs, method_es;

  std::cout << "Reading FCIdump file: " << filename<< "\n";

  // initialise hamiltonian
  method_gs = "CCSD";
  method_es = "EOM-CCSD";
  // init(filename, method_gs, hamiltonian); // regular electron only
  init_ccpol(filename, method_gs, hamiltonian); // add photon space
  std::cout << "HF energy: " << std::setprecision(12) 
            << energy_hf(hamiltonian.m2get(f_oo),hamiltonian.m4get(i_oooo))<< "\n";
// std::cout << "f_oo" << std::endl;
// hamiltonian.m2get(f_oo).print();
// std::cout << "i_oooo" << std::endl;
// hamiltonian.m4get(i_oooo).print();

#if 1 // CCSD
  // parse data
  if (method_gs == "CCSD" or method_gs == "CCS")
    ptampl->set(t1, container(hamiltonian.m2get(f_ov).get_space()));
  if (method_gs == "CCSD" or method_gs == "CCD")
    ptampl->set(t2, container(hamiltonian.m4get(i_oovv).get_space()));
  if (method_gs == "CCSD")
    problem.reset(new problem_ccsd(hamiltonian));
  if (method_gs == "CCD")
    problem.reset(new problem_ccd(hamiltonian));

  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<amplitudes<>>("DIIS");

  auto residual = *ptampl;
  // solver->set_verbosity(molpro::linalg::itsolv::Verbosity::None);
  // solver->set_convergence_threshold(1.0e-14);
  solver->solve(*ptampl, residual, *problem);
  solver->solution(*ptampl, residual);
  problem->energy(*ptampl);
  std::cout << method_gs << " correlation energy: " << std::setprecision(12) << problem->get_energy()<< "\n";
  // std::cout << method_gs << " total energy: " << std::setprecision(12) 
  //           << problem->get_energy() + energy_hf(hamiltonian.m2get(f_oo),hamiltonian.m4get(i_oooo))<< "\n";


  std::unique_ptr<amplitudes<>> ptampl2(new amplitudes());
  #if 0 // Excited State
  std::cout << "\n" << method_es<< "\n";
  size_t nroots(3);

  std::vector<amplitudes<>> v_rampl(nroots);
  std::unique_ptr<problem_eom> problem_es;
  
  problem_es.reset(new problem_eom_ccsd(hamiltonian, *ptampl));

  auto solver_es = molpro::linalg::itsolv::create_LinearEigensystem<amplitudes<>>("Davidson");
  auto residuals_es = v_rampl;
  solver_es->set_verbosity(molpro::linalg::itsolv::Verbosity::Iteration);
  solver_es->set_n_roots(nroots);
  solver_es->solve(v_rampl, residuals_es, *problem_es, true);
  problem_es->set_energy(solver_es->eigenvalues());

  #endif
  #endif

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "\nFinished computation at " << std::ctime(&end_time)
            << "Elapsed time: " << elapsed_seconds.count() << "s\n";
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<amplitudes<>, amplitudes<>>;
