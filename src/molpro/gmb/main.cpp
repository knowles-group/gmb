#include "init.h"
#include "problem_ccsd.h"
#include "problem_ccd.h"
#include "problem_eom-ccsd.h"
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
int main(int argc, char const *argv[])
{
  std::string filename{argv[0]};
  if (filename.find_last_of("/") != std::string::npos)
    filename.resize(filename.find_last_of("/"));
  filename += "/" + test_case + ".fcidump";
  auto start = std::chrono::system_clock::now();

  // std::unique_ptr<hamiltonian<>> hamiltonian;
  hamiltonian<> hamiltonian;
  amplitudes<> t_ampl;
  std::unique_ptr<problem_gen> problem;
  std::string method_gs, method_es;

  std::cout << "Reading FCIdump file: " << filename << std::endl;

  // initialise hamiltonian
  method_gs = "CCSD";
  method_es = "EOM-CCSD";
  init(filename, method_gs, hamiltonian);

  // parse data
  if (method_gs == "CCSD" or method_gs == "CCS")
    t_ampl.set(t1, hamiltonian.m2get(f_ov));
  if (method_gs == "CCSD" or method_gs == "CCD")
    t_ampl.set(t2, hamiltonian.m4get(i_oovv));
  if (method_gs == "CCSD")
    problem.reset(new problem_ccsd(hamiltonian));
  if (method_gs == "CCD")
    problem.reset(new problem_ccd(hamiltonian));

  auto solver = molpro::linalg::itsolv::create_NonLinearEquations<supercontainer<>>("DIIS");

  auto residual = t_ampl;
  solver->set_verbosity(molpro::linalg::itsolv::Verbosity::None);
  solver->solve(t_ampl, residual, *problem);
  solver->solution(t_ampl, residual);
  problem->energy(t_ampl);
  std::cout << method_gs << " energy: " << std::setprecision(8) << problem->get_energy() << std::endl;

  #if 1// Excited State
  std::cout << std::endl << method_es << std::endl;
  size_t nroots(2);
  supercontainer<> r_ampl;
  std::vector<supercontainer<>> v_rampl(nroots);
  std::vector<container<2>> vr1;
  std::unique_ptr<problem_eom> problem_es;
  for (size_t ir = 0; ir < nroots; ir++) {

    container r1(hamiltonian.m2get(f_ov).get_space());
    bbo::zero(r1);
    libtensor::block_tensor_wr_ctrl<2, double> ctrl(r1);
    libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
    for (libtensor::orbit_list<2, double>::iterator it = ol.begin();
         it != ol.end(); it++) {
      libtensor::index<2> bidx;
      ol.get_index(it, bidx);
      libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
      libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);
      const libtensor::dimensions<2> &tdims = blk.get_dims();
      double *ptr = tc.req_dataptr();
      if (ir == 1)
        if (bidx[0] != 1 or bidx[1] != 1) continue;
      ptr[0] = 1;
      tc.ret_dataptr(ptr);
      ctrl.ret_block(bidx);
      break;
    }
    v_rampl[ir].set("r1", r1);
    container r2(hamiltonian.m4get(i_oovv).get_space());
    bbo::zero(r2);
    v_rampl[ir].set("r2", r2);
  }
  
  problem_es.reset(new problem_eom_ccsd(hamiltonian, t_ampl));

  auto solver_es = molpro::linalg::itsolv::create_LinearEigensystem<supercontainer<>>("Davidson");
  auto residuals_es = v_rampl;
  solver_es->set_n_roots(nroots);
  solver_es->solve(v_rampl, residuals_es, *problem_es, false);
#endif

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "\nFinished computation at " << std::ctime(&end_time)
            << "Elapsed time: " << elapsed_seconds.count() << "s\n";
}

#include <molpro/linalg/itsolv/SolverFactory-implementation.h>
template class molpro::linalg::itsolv::SolverFactory<supercontainer<>, supercontainer<>>;
