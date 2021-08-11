#ifndef GMB_PROBLEM_EOM_H_
#define GMB_PROBLEM_EOM_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "amplitudes.h"
#include "hamiltonian.h"
#include <molpro/iostream.h>

class problem_eom : public molpro::linalg::itsolv::Problem<amplitudes<>> {
protected:
  std::vector<double> m_energy;    ///> energy
  mutable hamiltonian<> m_ham;     ///> Hamiltonian
  mutable amplitudes<> m_tampl;    ///> T amplitudes
  size_t m_nroots;    ///> energy
  mutable std::vector<std::unique_ptr<container<2>>> m_vr1;
public:
  using Problem::container_t;
  using Problem::value_t;
  problem_eom(const hamiltonian<> &ham, const amplitudes<> &tampl, const size_t &nroots)
  : m_ham(ham), m_tampl(tampl), m_nroots(nroots) {
    m_vr1.resize(nroots);
  }

  virtual ~problem_eom() {}

  void set_energy(std::vector<double> eigval) {m_energy = eigval;}

  std::vector<double> get_energy() const {return m_energy;}  
  
  virtual void create_guess(std::vector<amplitudes<>>& v_rampl) {};

  void character() const {

    for (size_t ir1 = 0; ir1 < m_vr1.size(); ir1++) {
      molpro::cout << "\nExcited state #" << ir1+1 
                << "\nExcitation energy: " << std::setprecision(5) << std::fixed << m_energy[ir1] << " Ha"
                << "\norbitals      transition type      amplitude\n";
      // molpro::cout << "r[" << ir1 << "]:\n";
      // m_vr1[ir1]->print();

      libtensor::block_tensor_rd_i<2, double> &bt(*m_vr1[ir1]);
      const libtensor::dimensions<2> &dims = bt.get_bis().get_dims();
      auto no = dims.get_dim(0);
      auto nv = dims.get_dim(1);

      libtensor::block_tensor_rd_ctrl<2, double> ctrl(*m_vr1[ir1]);

      libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
        size_t count{0};
      for (libtensor::orbit_list<2, double>::iterator it = ol.begin(); it != ol.end(); it++) {
        libtensor::index<2> bidx;
        ol.get_index(it, bidx);

        libtensor::dense_tensor_rd_i<2, double> &blk = ctrl.req_const_block(bidx);
        libtensor::dense_tensor_rd_ctrl<2, double> tc(blk);
        const libtensor::dimensions<2> &tdims = blk.get_dims();
        const double *ptr = tc.req_const_dataptr();
        size_t i{1}, a{1};
        for (size_t itdim = 0; itdim < tdims.get_size(); itdim++) {
          if (std::abs(ptr[itdim]) >  0.001) {
            i = 1+(count/nv);
            a = 1+itdim;
            molpro::cout << "o" << i  << " -> v" << a << "      ";
            for (size_t in = 0; in < 2; in++) {
              switch (bidx[in]) {
              case alpha: molpro::cout << "alpha ";
                break;
              case beta: molpro::cout << "beta  ";
                break;
              default: molpro::cout << "photon";
                break;
              }
              if (in == 0)
                molpro::cout << " -> ";
            }
            molpro::cout << "     " << std::setprecision(5) << std::fixed <<  ptr[itdim] << "\n";
          }
          ++count;
        }
        tc.ret_const_dataptr(ptr);
        ctrl.ret_const_block(bidx);
      }
    }
  };
};

#endif // GMB_PROBLEM_EOM_H_