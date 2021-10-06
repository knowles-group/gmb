#ifndef GMB_SRC_MOLPRO_GMB_PROBLEM_EOM_H_
#define GMB_SRC_MOLPRO_GMB_PROBLEM_EOM_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/iostream.h>
#include "amplitudes.h"
#include "hamiltonian.h"

class problem_eom : public molpro::linalg::itsolv::Problem<amplitudes<>> {

protected:
  std::vector<value_t> m_energy; ///> energy
  const std::shared_ptr<hamiltonian<>> m_ham;  ///> Hamiltonian
  const std::shared_ptr<amplitudes<>> m_tampl;  ///> T amplitudes

public:
  using Problem::container_t;
  using Problem::value_t;

  problem_eom(const std::shared_ptr<hamiltonian<>> &ham, const std::shared_ptr<amplitudes<>> &ptampl)
  : m_ham{ham}, m_tampl{ptampl} {}

  virtual ~problem_eom() = default;

  void set_energy(std::vector<value_t> eigval) {m_energy = eigval;}
  std::vector<value_t> get_energy() const {return m_energy;}  

  /**
   * @brief Give description of transitions in terms of orbitals
   * 
   */
  virtual void character(std::vector<container_t> &v_rampl) const {}

  /**
   * @brief Project out intruder triplet states
   * 
   * @param ampl r/l amplitudes
   */
  void singlet_projector(container_t &ampl) const {

    std::vector<size_t> v_no, v_nv;
    gmb::get_tensor_dimensions(ampl.m2get(r1), v_no, v_nv);
    std::vector<double> v_alpha, v_beta;
    
    // read r1  
    {
      constexpr size_t N = 2;
      libtensor::block_tensor_rd_ctrl<N, value_t> ctrl(ampl.m2get(r1));

      libtensor::orbit_list<N, value_t> ol(ctrl.req_const_symmetry());
      for (libtensor::orbit_list<N, value_t>::iterator it = ol.begin(); it != ol.end(); it++) {
        libtensor::index<N> bidx;
        ol.get_index(it, bidx);
        if (bidx[0] != bidx[1] || bidx[0] > beta)
            continue;
        libtensor::dense_tensor_rd_i<N, value_t> &blk = ctrl.req_const_block(bidx);
        libtensor::dense_tensor_rd_ctrl<N, value_t> tc(blk);
        const libtensor::dimensions<N> &tdims = blk.get_dims();
        const value_t *ptr = tc.req_const_dataptr();
        for (size_t offset = 0; offset < tdims.get_size(); offset++) {
          size_t i = 1+(offset/v_nv[bidx[1]]);
          size_t a = 1+offset-(offset/v_nv[bidx[1]])*v_nv[bidx[1]];
          switch (bidx[0]) {
          case alpha: v_alpha.push_back(ptr[offset]);
            break;
          case beta: v_beta.push_back(ptr[offset]);
            break;
          }
        }
        tc.ret_const_dataptr(ptr);
        ctrl.ret_const_block(bidx);
      }
        
      auto v_new = v_alpha;

      for (size_t i = 0; i < v_alpha.size(); i++) {
        v_new[i] = (v_alpha[i] + v_beta[i]) / 2;
      }

      // clean r1 
      {
        libtensor::block_tensor_wr_ctrl<N, value_t> ctrl(ampl.m2get(r1));
        libtensor::orbit_list<N, value_t> ol(ctrl.req_const_symmetry());
        for (libtensor::orbit_list<N, value_t>::iterator it = ol.begin(); it != ol.end(); it++) {
          libtensor::index<N> bidx;
          ol.get_index(it, bidx);
          if (bidx[0] != bidx[1] || bidx[0] > beta)
            continue;
          libtensor::dense_tensor_wr_i<N, value_t> &blk = ctrl.req_block(bidx);
          libtensor::dense_tensor_wr_ctrl<N, value_t> tc(blk);
          const libtensor::dimensions<N> &tdims = blk.get_dims();
          value_t *ptr = tc.req_dataptr();
          for (size_t offset = 0; offset < tdims.get_size(); offset++) {
              ptr[offset] = v_new[offset];
          }
          tc.ret_dataptr(ptr);
          ctrl.ret_block(bidx);
        }
      }
    }
  }


  };

#endif // GMB_PROBLEM_EOM_H_