#include "init.h"
#include "problem_ccsd.h"
#include "problem_ccd.h"
#include "hamiltonian.h"
#include "amplitudes.h"
#include <memory>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <chrono>
#include <ctime>

// test case
std::string test_case = "He-VDZ"; 
// std::string test_case = "H2O-VDZ"; 
// std::string test_case = "Li-VDZ-UHF"; 


using namespace bbo;
int main(int argc, char const *argv[]) {
  std::string filename{argv[0]};
  if (filename.find_last_of("/") != std::string::npos)
              filename.resize(filename.find_last_of("/"));
  filename += "/"+test_case+".fcidump";
  auto start = std::chrono::system_clock::now();

  // std::unique_ptr<hamiltonian<>> hamiltonian;
  hamiltonian<> hamiltonian;
  amplitudes<> t_ampl;
  std::unique_ptr<problem_gen> problem;
  std::string method;

  std::cout << "Reading FCIdump file: " << filename << std::endl;
  
  // initialise hamiltonian
  method = "CCSD";
  init(filename, method, hamiltonian);

  // parse data
  if (method == "CCSD" or method == "CCS")
    t_ampl.set(t1, hamiltonian.m2get(f_ov));
  if (method == "CCSD" or method == "CCD")
    t_ampl.set(t2, hamiltonian.m4get(i_oovv));
  if (method == "CCSD")
    problem.reset(new problem_ccsd(hamiltonian));
  if (method == "CCD")
    problem.reset(new problem_ccd(hamiltonian));


  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<supercontainer<>>("DIIS");

  auto residual = t_ampl;
  // solver->set_verbosity(molpro::linalg::itsolv::Verbosity::None);
  solver->solve(t_ampl, residual, *problem);
  solver->solution(t_ampl, residual);
  problem->energy(t_ampl);
  std::cout << method << " energy: " << std::setprecision(8) << problem->get_energy() << std::endl;

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "Finished computation at " << std::ctime(&end_time)
            << "Elapsed time: " << elapsed_seconds.count() << "s\n";
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<supercontainer<>, supercontainer<>>;
