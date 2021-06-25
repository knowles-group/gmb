#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "molpro/FCIdump.h"

#include "../src/molpro/gmb/init.cpp"
#include "../src/molpro/gmb/utils.cpp"
#include "../src/molpro/gmb/get_integral.cpp"
#include "../src/molpro/gmb/libtensor_utils.cpp"

#include "../src/molpro/gmb/problem_ccsd.h"
#include "../src/molpro/gmb/problem_eom-ccsd.h"
// #include "../src/molpro/gmb/hamiltonian.h"
#include "../src/molpro/gmb/amplitudes.h"

#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory.h>

// test case
std::string test_case = "He-VDZ";
// std::string test_case = "H2O-VDZ";
// std::string test_case = "Li-VDZ-UHF";
std::string filename(test_case+"/"+test_case+".fcidump");

std::unique_ptr<polariton> ppol;

int main(int argc, char* argv[]) {
  std::cout << "Test: " << test_case << std::endl;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


TEST(CCSD,energy) {
  hamiltonian<> hamiltonian;
  std::unique_ptr<amplitudes<>> ptampl(new amplitudes());
  std::unique_ptr<problem_gen> problem;
  std::string method_gs = "ccsd";

 // initialise hamiltonian
  gmb::init(filename, method_gs, hamiltonian);

  // initialise T amplitudes
  ptampl->set(t1, container(hamiltonian.m2get(f_ov).get_space()));
  ptampl->set(t2, container(hamiltonian.m4get(i_oovv).get_space()));
  problem.reset(new problem_ccsd(hamiltonian));

  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<amplitudes<>>("DIIS");

  auto residual = *ptampl;
  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::None);
  solver->set_convergence_threshold(1.0e-14);
  solver->solve(*ptampl, residual, *problem);
  solver->solution(*ptampl, residual);
  problem->energy(*ptampl);
  double energy = problem->get_energy();

  // get reference value and compare
  std::ifstream infile;
  std::string line;
  infile.open(test_case+"/energy_ccsd");
  getline(infile, line);
  double energy_ref = stod(line);
  ASSERT_NEAR(energy_ref, energy, 10E-13);
  infile.close();

#if 1
  // EOM-CCSD
  std::string method_es = "eom-ccsd";
  size_t nroots(3);

  std::vector<amplitudes<>> v_rampl(nroots);
  std::unique_ptr<problem_eom> problem_es;
  
  problem_es.reset(new problem_eom_ccsd(hamiltonian, *ptampl));

  auto solver_es = molpro::linalg::itsolv::create_LinearEigensystem<amplitudes<>>("Davidson");
  auto residuals_es = v_rampl;
  solver_es->set_verbosity(molpro::linalg::itsolv::Verbosity::None);
  solver_es->set_n_roots(nroots);
  solver_es->solve(v_rampl, residuals_es, *problem_es, true);
  problem_es->set_energy(solver_es->eigenvalues());
  energy = problem_es->get_energy()[1];

  // get reference value and compare
  infile.open(test_case+"/energy_eom-ccsd");
  getline(infile, line);
  energy_ref = stod(line);
  ASSERT_NEAR(energy_ref, energy, 10E-09);
  infile.close();
#endif
}



#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<amplitudes<>, amplitudes<>>;
