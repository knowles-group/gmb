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
#include "gmb.h"
#include <regex>

using namespace gmb;

extern molpro::Profiler prof;
extern "C" void general_many_body() { molpro::gmb::gmb();}

void molpro::gmb::gmb(const molpro::Options& options) {

  std::string filename = options.parameter("dump",std::string{""});
  auto expected_results = options.parameter("results",std::vector<double>{});
  std::vector<bool> found_expected_results(expected_results.size(),false);

  std::ios_base::sync_with_stdio(false);
  auto start = std::chrono::system_clock::now();

  // parse arguments
  auto method = options.parameter("method","eom-ccsd");
  auto nroots = options.parameter("states", 3);
  auto ncav = options.parameter("polariton_modes", 0);

  std::vector<std::shared_ptr<polariton>> v_ppol(ncav);
  if (ncav > 0) {
    std::vector<int> v_option_polariton_nmax(ncav, 1); 
    v_option_polariton_nmax = options.parameter("polariton_nmax", v_option_polariton_nmax);
    std::vector<double> v_option_polariton_gamma(ncav, 0.01); 
    v_option_polariton_gamma = options.parameter("polariton_gamma", v_option_polariton_gamma);
    std::vector<double> v_option_polariton_omega(ncav, 1.028);
    v_option_polariton_omega = options.parameter("polariton_omega", v_option_polariton_omega);
    for (size_t i = 0; i < ncav; i++) {
      v_ppol[i] = std::make_shared<polariton>(v_option_polariton_nmax[i],
                                              v_option_polariton_gamma[i],
                                              v_option_polariton_omega[i]);
      v_ppol[i]->fname_dm = options.parameter(
          "dipole", std::regex_replace(
                        filename, std::regex{"\\.[_[:alnum:]]*$"}, ".dm"));
      v_ppol[i]->fname_sm = std::regex_replace(
                        filename, std::regex{"\\.[_[:alnum:]]*$"}, ".sm");
    }
  }

  // print arguments
  std::cout << "Required calculation: " << "\n\n";
  check_file(filename, "fcidump");
  std::cout << " fcidump = " << filename << "\n";
  std::cout << " method = " << method << "\n";
  std::cout << " roots = " << nroots << "\n";
  if (v_ppol.size() > 0) {
    check_file(v_ppol[0]->fname_dm, "dipole");
    std::cout << " dipole file = " << v_ppol[0]->fname_dm << "\n";
    check_file(v_ppol[0]->fname_sm, "second moment of charges");
    std::cout << " second moment of charges file = " << v_ppol[0]->fname_sm << "\n";

    std::cout << "\nPolariton parameters: \n";
    std::cout << "modes: " << ncav << "\n";
    for (size_t i = 0; i < ncav; i++) {
      std::cout << "\n mode " << i << "\n";
      std::cout << " nmax = " << v_ppol[i]->nmax << "\n";
      std::cout << " gamma = " << v_ppol[i]->gamma << "\n";
      std::cout << " omega = " << v_ppol[i]->omega << "\n\n";
    }
  }

  hamiltonian<> ham;
  std::unique_ptr<amplitudes<>> ptampl{std::make_unique<amplitudes<>>()};
  std::unique_ptr<problem_gen> problem;
  std::string method_gs, method_es;

  // initialise hamiltonian
  init(filename, method, ham, v_ppol);

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
    problem = std::make_unique<problem_ccsd>(ham);
  if (method.find("ccd") != std::string::npos)
    problem = std::make_unique<problem_ccd>(ham);

  // set solver
  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<amplitudes<>>("DIIS", "max_size_qspace=10");
  auto residual = *ptampl;

    // solver options
  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Iteration);
  // solver->set_max_iter(110);
  // solver->set_convergence_threshold(1.0e-7);
  solver->solve(*ptampl, residual, *problem);
  solver->solution(*ptampl, residual);
  problem->energy(*ptampl);

  // print results
  std::cout << *problem << " correlation energy: " << std::setprecision(12) << problem->get_energy()<< "\n";
  std::cout << *problem  << " total energy: " << std::setprecision(13)
            << problem->get_energy() + hf_energy<< "\n";
  for (int i=0; i<expected_results.size(); ++i)
    if (std::abs(problem->get_energy()+hf_energy-expected_results[i])<1e-9) found_expected_results[i]=true;


#if 1 // Excited State

  // set EOM-CCSD amplitudes
  std::unique_ptr<amplitudes<>> prampl{std::make_unique<amplitudes<>>()};
  prampl->set(r1, container(ptampl->m2get(t1).get_space()));
  prampl->set(r2, container(ptampl->m4get(t2).get_space()));
  std::vector<amplitudes<>> v_rampl(nroots, *prampl);
  std::unique_ptr<problem_eom> problem_es;
  problem_es = std::make_unique<problem_eom_ccsd>(ham, *ptampl);
  
  std::cout << "\n" << *problem_es << "\n";

  // set solver 
  auto solver_es = molpro::linalg::itsolv::create_LinearEigensystem<amplitudes<>>("Davidson");
  auto residuals_es = v_rampl;
  
  // set options
  solver_es->set_verbosity(molpro::linalg::itsolv::Verbosity::Iteration);
  solver_es->set_n_roots(nroots);
  solver_es->set_convergence_threshold(1.0e-7);

  // solve
  solver_es->solve(v_rampl, residuals_es, *problem_es, true);
  problem_es->set_energy(solver_es->eigenvalues());
  for (const auto& ev : solver_es->eigenvalues())
    for (int i=0; i<expected_results.size(); ++i)
      if (std::abs(ev-expected_results[i])<1e-9) found_expected_results[i]=true;
  auto energies = problem_es->get_energy();

  // print results
  std::cout << "\n" << *problem_es << " excitation energies (Ha) \n";
  for (auto &i : energies)
    std::cout << i << " \n";

  #endif
  #endif

  // std::cout << *pprof << std::endl;
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "\nFinished computation at " << std::ctime(&end_time)
            << "Elapsed time: " << elapsed_seconds.count() << "s\n";
  for (int i=0; i<expected_results.size(); ++i)
    if (not found_expected_results[i])
      throw std::runtime_error("Did not match expected result "+std::to_string(expected_results[i]));

}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
// template class molpro::linalg::itsolv::SolverFactory<amplitudes<>, amplitudes<>>;
template class molpro::linalg::itsolv::SolverFactory<amplitudes<>, amplitudes<>>;
