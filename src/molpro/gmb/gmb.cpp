#include <molpro/Options.h>
#include <molpro/iostream.h>
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory.h>

#include <chrono>
#include <ctime>
#include <memory>
#include <regex>

#include "amplitudes.h"
#include "constants.h"
#include "expressions/energy_hf.h"
#include "gmb.h"
#include "hamiltonian.h"
#include "init.h"
#include "run_method.h"

using namespace gmb;

extern molpro::Profiler prof;
extern "C" void general_many_body(int64_t &nstate, double *energies) {
  auto ev = molpro::gmb::gmb();
  nstate = ev.size();
  for (size_t i = 0; i < nstate; ++i)
    energies[i] = ev[i];
}

std::vector<double> molpro::gmb::gmb(const molpro::Options &options) {

  std::vector<double> all_energies;
  #if 1
  const auto filename = options.parameter("dump", std::string{""});
  const auto expected_results = options.parameter("results", std::vector<double>{});
  std::vector<bool> found_expected_results(expected_results.size(), false);

  std::ios_base::sync_with_stdio(false);
  const auto start = std::chrono::system_clock::now();

  // parse arguments
  const auto method = options.parameter("method", "eom-ccsd");
  const auto nroots = options.parameter("states", 3);
  const auto es_conv = options.parameter("es_conv", 1e-5);
  const auto ncav = options.parameter("polariton_modes", 0);
  const auto nvib = options.parameter("vibration_modes", 0);

  std::vector<std::shared_ptr<polariton>> v_ppol(ncav);
  if (ncav > 0) {
    const auto self_energy = options.parameter("self_energy", true);
    const auto coupling = options.parameter("coupling", true);
    std::vector<int> v_option_polariton_nmax(ncav, 1);
    v_option_polariton_nmax =
        options.parameter("polariton_nmax", v_option_polariton_nmax);
    std::vector<double> v_option_polariton_gamma(ncav, 0.01);
    v_option_polariton_gamma =
        options.parameter("polariton_gamma", v_option_polariton_gamma);
    std::vector<double> v_option_polariton_omega(ncav, 1.0);
    v_option_polariton_omega =
        options.parameter("polariton_omega", v_option_polariton_omega);
    for (size_t i = 0; i < ncav; i++) {
      v_ppol[i] = std::make_shared<polariton>(v_option_polariton_nmax[i],
                                              v_option_polariton_gamma[i],
                                              v_option_polariton_omega[i],
                                              self_energy,
                                              coupling);
      v_ppol[i]->fname_dm = options.parameter(
          "dipole",
          std::regex_replace(filename, std::regex{"\\.[_[:alnum:]]*$"}, ".dm"));
      v_ppol[i]->fname_sm =
          std::regex_replace(filename, std::regex{"\\.[_[:alnum:]]*$"}, ".sm");
    }
  }

  std::vector<std::unique_ptr<vibration>> v_pvib(nvib);

  if (nvib > 0) {
    std::vector<int> v_option_vibration_nmax(nvib, 1);
    v_option_vibration_nmax =
        options.parameter("vibration_nmax", v_option_vibration_nmax);
    std::vector<double> v_option_vibration_omega(nvib, 0.01);
    v_option_vibration_omega =
        options.parameter("vibration_omega", v_option_vibration_omega);
    std::vector<double> v_option_vibration_damping(nvib, 0.0);
    v_option_vibration_damping =
        options.parameter("vibration_damping", v_option_vibration_damping);
    for (size_t i = 0; i < nvib; i++) {
      v_pvib[i] = std::make_unique<vibration>(v_option_vibration_nmax[i],
                                              v_option_vibration_omega[i]);
      v_pvib[i]->integral_files[0] = options.parameter(
          "c",
          std::regex_replace(filename, std::regex{"\\.[_[:alnum:]]*$"}, ".c"));
      v_pvib[i]->integral_files[1] =
          std::regex_replace(filename, std::regex{"\\.[_[:alnum:]]*$"}, ".a");
      v_pvib[i]->integral_files[2] =
          std::regex_replace(filename, std::regex{"\\.[_[:alnum:]]*$"}, ".pi");
      v_pvib[i]->integral_files[3] =
          std::regex_replace(filename, std::regex{"\\.[_[:alnum:]]*$"}, ".PI");
    }
  }

  
  
  // print arguments
  molpro::cout << "\nRequired calculation: "
               << "\n\n";
  check_file(filename, "fcidump");
  molpro::cout << " fcidump = " << filename << "\n";
  molpro::cout << " method = " << method << "\n";
  molpro::cout << " roots = " << nroots << "\n\n";

  if (!v_ppol.empty() ) {

    molpro::cout << "\nPolariton parameters: \n";

    check_file(v_ppol[0]->fname_dm, "dipole");
    molpro::cout << " dipole file = " << v_ppol[0]->fname_dm << "\n";
    check_file(v_ppol[0]->fname_sm, "second moment of charges");
    molpro::cout << " second moment of charges file = " << v_ppol[0]->fname_sm
                 << "\n";
    
    molpro::cout << "nº modes: " << ncav << "\n";
    for (size_t i = 0; i < ncav; i++) {
      molpro::cout << "\n mode " << i+1 << "\n";
      molpro::cout << "  nmax = " << v_ppol[i]->nmax << "\n";
      molpro::cout << "  gamma = " << v_ppol[i]->gamma << "\n";
      molpro::cout << "  omega = " << v_ppol[i]->omega << "\n\n";
    }
  }

  if (!v_pvib.empty() ) {
    molpro::cout << "\nVibrational parameters: \n";

    check_file(v_pvib[0]->integral_files[0], "c");
    check_file(v_pvib[0]->integral_files[1], "a");
    check_file(v_pvib[0]->integral_files[2], "pi");
    check_file(v_pvib[0]->integral_files[3], "PI");

    molpro::cout << "\nIntegral files:\n"
              << " constant c file is : " << v_pvib[0]->integral_files[0] << "\n"
              << " A matrix file is : " << v_pvib[0]->integral_files[1] << "\n"
              << " pi matrix file is : " << v_pvib[0]->integral_files[2] << "\n"
              << " PI matrix file is : " << v_pvib[0]->integral_files[3] << "\n";

    molpro::cout << "\nnº modes: " << nvib << "\n";
    for (size_t i = 0; i < nvib; i++) {
      molpro::cout << "\n mode " << i+1 << ":\n"
                   << "  nmax = " << v_pvib[i]->nmax << "\n"
                   << "  omega = " << v_pvib[i]->omega << "\n\n";
    }
  }

  // initialise hamiltonian
  hamiltonian<> ham;
  init(filename, method, ham, v_ppol, v_pvib);

#if 1 // GS
  const auto vnn = get_integral(filename);
  auto hf_energy = vnn + energy_hf(ham.m2get(f_oo),ham.m4get(i_oooo));

  #if 1 // self-energy
  for (size_t i = 0; i < ncav; i++) {
    if (v_ppol[i]->self_energy) {
      const auto rnuc = get_integral(v_ppol[i]->fname_dm);
      molpro::cout << " rnuc = " << rnuc << "\n";
      hf_energy += v_ppol[i]->gamma*v_ppol[i]->gamma*v_ppol[i]->omega*rnuc*rnuc;
    }
  }
  #endif
  molpro::cout << "\nHF energy: " << std::setprecision(12) << hf_energy << "\n\n";

  
  if (!(method.find("hf") != std::string::npos)) {

    std::unique_ptr<problem_gen> problem;
    std::unique_ptr<amplitudes<>> ptampl{std::make_unique<amplitudes<>>()};
    
    run_gs(ham, method, problem, ptampl);
    
    // print results
    molpro::cout << *problem << " correlation energy: " << std::setprecision(12) << problem->get_energy()<< "\n";
    const auto ccsd_energy = problem->get_energy() + hf_energy;
    all_energies.push_back(ccsd_energy);
    molpro::cout << *problem  << " total energy: " << std::setprecision(13)
                 << ccsd_energy << "\n";
    for (int i=0; i<expected_results.size(); ++i)
      if (std::abs(problem->get_energy()+hf_energy-expected_results[i])<1e-10) found_expected_results[i]=true;

#if 1 // ES
    if (method.find("eom") != std::string::npos) {

      std::unique_ptr<problem_eom> problem_es;
      // problem_es = std::make_unique<problem_eom_ccsd>(ham, *ptampl, v_ppol);
      run_eom(ham, method, problem_es, ptampl, nroots, es_conv, ccsd_energy, v_ppol);

      for (const auto& ev : problem_es->get_energy())
        for (int i=0; i<expected_results.size(); ++i)
          if (std::abs(ev-expected_results[i])<1e-10) found_expected_results[i]=true;
      auto energies = problem_es->get_energy();

      // print results
      molpro::cout << "\n          Excitation energy                   Total energy  \n";
      molpro::cout << "        (Ha)            (eV)              (Ha)            (eV)  \n";
      for (auto &i : energies) {
        all_energies.push_back(ccsd_energy+i);
        molpro::cout << std::fixed << std::setw(14) << std::setprecision(7) << i << "  "
                     << std::fixed << std::setw(14) << std::setprecision(6) << i*inverse_electron_volt << "      "
                     << std::fixed << std::setw(14) << std::setprecision(10) << ccsd_energy+i << "  "
                     << std::fixed << std::setw(14) << std::setprecision(6) << (ccsd_energy+i)*inverse_electron_volt 
                     << "\n";
      }
      molpro::cout << std::endl;
    }
  #endif
  }
  #endif

  const auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  molpro::cout << "\nFinished computation at " << std::ctime(&end_time)
               << "Elapsed time: " << elapsed_seconds.count() << "s\n\n";
  for (int i = 0; i < expected_results.size(); ++i)
    if (not found_expected_results[i])
      throw std::runtime_error("Did not match expected result " +
          std::to_string(expected_results[i]));
  #endif 
  return all_energies;
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<amplitudes<>, amplitudes<>>;
