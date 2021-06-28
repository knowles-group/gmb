#include "init.h"
#include "problem_ccsd.h"
#include "problem_ccd.h"
#include "problem_eom-ccsd.h"
#include "hamiltonian.h"
#include "amplitudes.h"
#include "expressions/energy_hf.h"
#include "utils.h"
#include <memory>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory.h>
#include <chrono>
#include <ctime>
#include <molpro/Options.h>
#include "expressions/diag_ov.h"

using namespace gmb;

// test case
std::string filename;
std::string test_case = "hubbard"; 
std::unique_ptr<polariton> ppol;

int main(int argc, char *argv[]) {
  molpro::Options options("gmb", argc, argv);
  filename = options.parameter("dump",std::string{""});
  std::ios_base::sync_with_stdio(false);
//  filename = argv[0];
//  if (filename.find_last_of("/") != std::string::npos)
//    filename.resize(filename.find_last_of("/"));
//  filename += "/"+test_case+"/"+test_case+".fcidump";
  auto start = std::chrono::system_clock::now();

  // parse arguments
  polariton pol;
  std::string dump{filename}, method{"eom-ccsd"};
  // syms_t states(8);
  int nroots{7};
  for (int i = 0; i < argc; ++i) {
    std::string arg{argv[i]};
//    if (arg.substr(0,arg.find('=')) == "dump") filename = arg.substr(arg.find('=')+1);
    if (arg.substr(0,arg.find('=')) == "method") method = arg.substr(arg.find('=')+1);
    if (arg.substr(0,arg.find('=')) == "states") {
      std::stringstream ss;
      ss << arg.substr(arg.find('=')+1);
      ss >> nroots;
    }
    if (arg.substr(0,arg.find('=')) == "polariton") {
      std::stringstream ss{arg.substr(arg.find('[')+1, arg.find(']')-1-arg.find('['))};
      polariton pol;
      get_polval(ss,pol.nmax);
      get_polval(ss,pol.gamma);
      get_polval(ss,pol.omega);
      ppol = std::make_unique<polariton>(pol.nmax,pol.gamma,pol.omega);
    }
  }
  // ppol = std::make_unique<polariton>(1,0.01,1.028);

  std::cout << "Required calculation: " << "\n";
  std::cout << "dump = " << filename << "\n";
  std::cout << "method = " << method << "\n";
  std::cout << "roots = " << nroots << "\n";
  if (ppol != nullptr) {
    std::cout << "polariton parameters:" << std::endl;
    std::cout << "nmax = " << ppol->nmax << std::endl;
    std::cout << "gamma = " << ppol->gamma << std::endl;
    std::cout << "omega = " << ppol->omega << std::endl;
  }

  hamiltonian<> ham;
  std::unique_ptr<amplitudes<>> ptampl(std::make_unique<amplitudes<>>());
  std::unique_ptr<problem_gen> problem;
  std::string method_gs, method_es;

  // initialise hamiltonian
  init(filename, method, ham); // add photon space
    

#if 1 // CCSD
  auto vnn = get_integral(filename);
  auto hf_energy = vnn + energy_hf(ham.m2get(f_oo),ham.m4get(i_oooo));
  std::cout << "HF energy: " << std::setprecision(12) << hf_energy << "\n";


  // set CCSD amplitudes
  if (method.find("ccs") != std::string::npos)
    ptampl->set(t1, container(ham.m2get(f_ov).get_space()));
  if ((method.find("ccsd") != std::string::npos) || (method.find("ccd") != std::string::npos))
    ptampl->set(t2, container(ham.m4get(i_oovv).get_space()));

  // set problem
  if (method.find("ccsd") != std::string::npos)
    problem.reset(new problem_ccsd(ham));
  if (method.find("ccd") != std::string::npos)
    problem.reset(new problem_ccd(ham));

  // set solver
  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<amplitudes<>>("DIIS");
  auto residual = *ptampl;

  // solver options
  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Iteration);
  // solver->set_max_iter(110);
  // solver->set_convergence_threshold(1.0e-14);
  solver->solve(*ptampl, residual, *problem);
  solver->solution(*ptampl, residual);
  problem->energy(*ptampl);

  // print results
  std::cout << *problem << " correlation energy: " << std::setprecision(12) << problem->get_energy()<< "\n";
  std::cout << *problem  << " total energy: " << std::setprecision(13) 
            << problem->get_energy() + hf_energy<< "\n";
  
  #if 1 // Excited State
  std::vector<amplitudes<>> v_rampl(nroots);
  std::unique_ptr<problem_eom> problem_es;
  
  problem_es.reset(new problem_eom_ccsd(ham, *ptampl));
  std::cout << "\n" << *problem_es << "\n";
  auto solver_es = molpro::linalg::itsolv::create_LinearEigensystem<amplitudes<>>("Davidson");
  auto residuals_es = v_rampl;
  solver_es->set_verbosity(molpro::linalg::itsolv::Verbosity::Iteration);
  solver_es->set_n_roots(nroots);
  solver_es->solve(v_rampl, residuals_es, *problem_es, true);
  problem_es->set_energy(solver_es->eigenvalues());
  auto energies = problem_es->get_energy();
  std::cout << "\n" << *problem_es << " excitation energies (Ha) \n";  
  for (auto &i : energies)
    std::cout << i << " \n";  

  #endif
  #endif

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "\nFinished computation at " << std::ctime(&end_time)
            << "Elapsed time: " << elapsed_seconds.count() << "s\n";
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
// template class molpro::linalg::itsolv::SolverFactory<amplitudes<>, amplitudes<>>;
template class molpro::linalg::itsolv::SolverFactory<amplitudes<>, amplitudes<>>;
