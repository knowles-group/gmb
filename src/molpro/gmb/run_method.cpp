#include "run_method.h"
#include "problem_eom-ccsd.h"
#include "problem_eom-ccsd-left.h"

#include "problem_ccsd.h"
#include "problem_ccd.h"
#include "problem_mp2.h"

#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/linalg/itsolv/SolverFactory.h>

#include <molpro/iostream.h>

void run_gs(hamiltonian<> &ham, const std::string &method, std::unique_ptr<problem_gen> &problem, std::shared_ptr<amplitudes<>> &ptampl) {

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

void run_eom(const std::shared_ptr<hamiltonian<>> &pham, const std::string &method, 
             std::unique_ptr<problem_eom> &problem, const std::shared_ptr<amplitudes<>> &ptampl, 
             const size_t &nroots, const double& es_conv) {
  molpro::cout << "\nRunning EOM-CCSD" << std::endl;
  
  supercontainer<> intermediates;
  auto pint = std::make_shared<supercontainer<>>(intermediates);
  eom_intermediates(pham, ptampl, pint);

  // set EOM-CCSD amplitudes
  problem = std::make_unique<problem_eom_ccsd>(pham, ptampl, pint);

  auto prampl = std::make_unique<amplitudes<>>();
  prampl->set(r1, container(ptampl->m2get(t1).get_space()));
  prampl->set(r2, container(ptampl->m4get(t2).get_space()));
  std::vector<amplitudes<>> v_rampl(nroots, *prampl);
  
  // set solver 
  {
    auto solver = molpro::linalg::itsolv::create_LinearEigensystem<amplitudes<>>("Davidson");
    auto residuals_es = v_rampl;

    // set options
    solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Iteration);
    solver->set_n_roots(nroots);
    solver->set_convergence_threshold(es_conv);

    // solve
    solver->solve(v_rampl, residuals_es, *problem, true);
    problem->set_energy(solver->eigenvalues());
    std::vector<int> v_nroots(nroots);
    for (size_t i = 0; i < nroots; i++)
      v_nroots[i] = i;
    solver->solution(v_nroots, v_rampl, residuals_es);
    problem->character(v_rampl);
  }
  auto energies = problem->get_energy();

  bool left{true};
  if (left) {
    // set EOM-CCSD left amplitudes
    molpro::cout << "\nSolving left EOM-CCSD\n";
    problem = std::make_unique<problem_eom_ccsd_left>(pham, ptampl, pint);

    auto v_lampl = v_rampl;

    // set solver 
    auto solver = molpro::linalg::itsolv::create_LinearEigensystem<amplitudes<>>("Davidson");
    auto residuals = v_lampl;

    // set options
    solver->set_verbosity(molpro::linalg::itsolv::Verbosity::Iteration);
    solver->set_n_roots(nroots);
    solver->set_convergence_threshold(es_conv);

    // solve
    solver->solve(v_lampl, residuals, *problem, false);
    // problem->set_energy(solver->eigenvalues());
    // std::vector<int> v_nroots(nroots);
    // for (size_t i = 0; i < nroots; i++)
    //   v_nroots[i] = i;
    // solver->solution(v_nroots, v_lampl, residuals);
    // problem->character(v_lampl);

    // auto energies = problem->get_energy();
  }
}

  void eom_intermediates(const std::shared_ptr<hamiltonian<>> &pham, 
                         const std::shared_ptr<amplitudes<>> &tampl, 
                         std::shared_ptr<supercontainer<>> &intermediates) {

    auto tau = ccsd_tau(tampl->m2get(t1), tampl->m4get(t2));
    auto if_oo = eom_ccsd_if_oo(tampl->m2get(t1), tampl->m4get(t2), pham->m2get(f_oo), pham->m2get(f_ov),
                   pham->m4get(i_ooov), pham->m4get(i_oovv));    
    auto if_ov = eom_ccsd_if_ov(tampl->m2get(t1), pham->m2get(f_ov), pham->m4get(i_oovv)); 
    auto if_vv = eom_ccsd_if_vv(tampl->m2get(t1), tampl->m4get(t2), pham->m2get(f_ov), pham->m2get(f_vv), 
                   pham->m4get(i_oovv), pham->m4get(i_ovvv));

    auto iw_oooo = eom_ccsd_iw_oooo(tampl->m2get(t1), tau, 
                    pham->m4get(i_oooo), pham->m4get(i_ooov), pham->m4get(i_oovv));    
    auto iw_ooov = eom_ccsd_iw_ooov(tampl->m2get(t1), tampl->m4get(t2), tau, if_ov,
                    iw_oooo, pham->m4get(i_ooov), pham->m4get(i_oovv), pham->m4get(i_ovov), pham->m4get(i_ovvv));
    auto iw_ovov = eom_ccsd_iw_ovov(tampl->m2get(t1), tampl->m4get(t2), 
                    pham->m4get(i_ooov), pham->m4get(i_oovv), pham->m4get(i_ovov), pham->m4get(i_ovvv));
    auto iw_vvvv = eom_ccsd_iw_vvvv(tampl->m2get(t1), tau, pham->m4get(i_oovv), pham->m4get(i_ovvv), pham->m4get(i_vvvv));
    auto iw_ovvv = eom_ccsd_iw_ovvv(tampl->m2get(t1), tampl->m4get(t2), tau, if_ov,  
                    pham->m4get(i_ooov), pham->m4get(i_oovv), pham->m4get(i_ovov), pham->m4get(i_ovvv), iw_vvvv);
    auto iw2_ooov = eom_ccsd_iw2_ooov(tampl->m2get(t1), pham->m4get(i_ooov), pham->m4get(i_oovv));
    auto iw2_ovvv = eom_ccsd_iw2_ovvv(tampl->m2get(t1), pham->m4get(i_oovv), pham->m4get(i_ovvv));
 
    intermediates->set("tau", tau);
    intermediates->set("if_oo", if_oo);
    intermediates->set("if_ov", if_ov);
    intermediates->set("if_vv", if_vv);
    intermediates->set("iw_oooo", iw_oooo);
    intermediates->set("iw_ooov", iw_ooov);
    intermediates->set("iw_ovov", iw_ovov);
    intermediates->set("iw_ovvv", iw_ovvv);
    intermediates->set("iw_vvvv", iw_vvvv);
    intermediates->set("iw2_ooov", iw2_ooov);
    intermediates->set("iw2_ovvv", iw2_ovvv);   
  }
  
#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<amplitudes<>, amplitudes<>>;