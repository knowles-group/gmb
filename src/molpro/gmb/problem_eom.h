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
    constexpr size_t N =2;
    constexpr double inverse_electron_volt{27.211386245988};
    for (size_t ir1 = 0; ir1 < m_vr1.size(); ir1++) {
      // molpro::cout << "r[" << ir1 << "]:\n";
      // m_vr1[ir1]->print();
      molpro::cout << "\nExcited state #" << ir1+1 
                << "\nExcitation energy = " << std::setprecision(5) << std::fixed 
                << m_energy[ir1] << " Ha  = "
                << m_energy[ir1]*inverse_electron_volt << " eV"
                << "\nocc -> vir     amplitude\n";

      libtensor::block_tensor_rd_i<2, double> &bt(*m_vr1[ir1]);

      // total dimensions
      const libtensor::dimensions<2> &dims = bt.get_bis().get_dims();
      auto no = dims.get_dim(0);
      auto nv = dims.get_dim(1);

      auto bis = bt.get_bis();

      size_t maxtyp = 0;
      for(size_t i = 0; i < N; i++) {
          auto typ = bis.get_type(i);
          if(typ > maxtyp) maxtyp = typ;
      }
      std::vector<size_t> v_no;
      std::vector<size_t> v_nv;

    // occupied
    const libtensor::split_points &spl_o = bis.get_splits(0);
    for (size_t i = 0; i < spl_o.get_num_points(); i++){
      if (i == 0)
        v_no.push_back(spl_o[i]);
      else 
        v_no.push_back(spl_o[i]-spl_o[i-1]);
    }
    v_no.emplace_back(no-std::accumulate(v_no.cbegin(),v_no.cend(),0));

    // virtual
    const libtensor::split_points &spl_v = bis.get_splits(1);
    for (size_t i = 0; i < spl_v.get_num_points(); i++) {
      if (i == 0)
        v_nv.push_back(spl_v[i]);
      else 
        v_nv.push_back(spl_v[i]-spl_v[i-1]);
    }
    v_nv.emplace_back(nv-std::accumulate(v_nv.cbegin(),v_nv.cend(),0));    

    std::vector<std::vector<size_t>> n_ne{v_no,v_nv};
      libtensor::block_tensor_rd_ctrl<2, double> ctrl(*m_vr1[ir1]);

      libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
        size_t count{0};
      for (libtensor::orbit_list<2, double>::iterator it = ol.begin(); it != ol.end(); it++) {
        libtensor::index<2> bidx;
        ol.get_index(it, bidx);
        const libtensor::dimensions<2> &bd = bt.get_bis().get_block_dims(bidx);
        libtensor::dense_tensor_rd_i<2, double> &blk = ctrl.req_const_block(bidx);
        libtensor::dense_tensor_rd_ctrl<2, double> tc(blk);
        const libtensor::dimensions<2> &tdims = blk.get_dims();
        const double *ptr = tc.req_const_dataptr();
        for (size_t offset = 0; offset < tdims.get_size(); offset++) {
          if (std::abs(ptr[offset]) >  0.001) {
            size_t i = 1+(offset/v_nv[bidx[1]]);
            molpro::cout << "o" << i;
            switch (bidx[0]) {
            case alpha: molpro::cout << "a";
              break;
            case beta: molpro::cout << "b";
              break;
            default: molpro::cout << "p";
              break;
            }
            size_t a = 1+offset-(offset/v_nv[bidx[1]])*v_nv[bidx[1]];
            molpro::cout << " -> v" << a;
            switch (bidx[1]) {
            case alpha: molpro::cout << "a";
              break;
            case beta: molpro::cout << "b";
              break;
            default: molpro::cout << "p";
              break;
            }
            molpro::cout << "     " << std::setprecision(5) << std::fixed <<  ptr[offset] << "\n";
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