#ifndef GMB_SRC_MOLPRO_GMB_PROBLEM_EOM_CCSD_H_
#define GMB_SRC_MOLPRO_GMB_PROBLEM_EOM_CCSD_H_
#include "constants.h"
#include "expressions/eom-ccsd/eom-ccsd.h"
#include "expressions/update.h"
#include "expressions/ccsd/energy.h"
#include "problem_eom.h"
#include "utils.h"

#include <chrono>
#include <ctime>
/**
 * @brief EOM-CCSD
 * 
 */

class problem_eom_ccsd : public problem_eom {
private:
  const std::shared_ptr<supercontainer<>> m_int;       ///> intermediates
public:
  problem_eom_ccsd(const std::shared_ptr<hamiltonian<>> &ham, 
                   const std::shared_ptr<amplitudes<>> &ptampl,
                   const std::shared_ptr<supercontainer<>> &intermediates)
  : problem_eom{ham, ptampl}, m_int(intermediates) {}


  bool diagonals(container_t &d) const override {
    auto d_ov = diag_ov(m_int->m2get("if_oo"), m_int->m2get("if_ov"), m_int->m2get("if_vv"));
    auto d_oovv = diag_oovv(d_ov, m_ham->m4get(i_oovv));
    auto id_oo = id_xx(m_ham->m2get(f_oo));
    auto id_vv = id_xx(m_ham->m2get(f_vv));
    auto d1 = hbar_ov(d_ov);
    auto d2 = hbar_oovv(id_oo, id_vv, d_oovv);
    d.set(r1,d1);
    d.set(r2,d2);
    return true;
  }

  void precondition(const VecRef<container_t>& residual, 
                    const std::vector<value_t>& shift, 
                    const container_t& diagonals) const override {
    for (int k = 0; k < residual.size(); k++) {
      auto &a = residual[k].get();   
      auto &d = const_cast<container_t&> (diagonals);  
      auto a1 = precond_ov(a.m2get(r1), d.m2get(r1), shift[k]);
      auto id_oo = id_xx(m_ham->m2get(f_oo));
      auto id_vv = id_xx(m_ham->m2get(f_vv));
      auto a2 = precond_oovv(a.m4get(r2),d.m4get(r2), shift[k]);
      a.set(r1,a1);
      a.set(r2,a2);
    }
  }

  void action(const CVecRef<container_t> &parameters, const VecRef<container_t> &actions) const override {
  
    for (int k = 0; k < parameters.size(); k++) {

      auto &ccp = const_cast<container_t&> (parameters[k].get());     
      auto &a = actions[k].get();  

      singlet_projector(ccp);

      // compute intermediates
      auto ir1_vv = eom_ccsd_ir1_vv(m_ham->m2get(f_vv),ccp.m2get(r1),m_int->m4get("iw2_ovvv"));          
      auto ir1_oo = eom_ccsd_ir1_oo(m_ham->m2get(f_oo),ccp.m2get(r1),m_int->m4get("iw2_ooov"));          
      auto ir2_oo = eom_ccsd_ir2_oo(m_ham->m2get(f_oo),ccp.m4get(r2),m_ham->m4get(i_oovv));          
      auto ir2_vv = eom_ccsd_ir2_vv(m_ham->m2get(f_vv),ccp.m4get(r2),m_ham->m4get(i_oovv));          

      // compute r1
      {
        auto r1_new = eom_ccsd_r1(ccp.m2get(r1), ccp.m4get(r2), m_int->m2get("if_oo"), m_int->m2get("if_ov"), m_int->m2get("if_vv"),  
                      m_int->m4get("iw_ovov"), m_int->m4get("iw2_ooov"), m_int->m4get("iw2_ovvv"));
        a.set(r1, r1_new);
      }
      // compute r2
      {
        auto r2_new = eom_ccsd_r2(ccp.m2get(r1), ccp.m4get(r2), m_tampl->m4get(t2), m_int->m2get("if_oo"), m_int->m2get("if_vv"),  
                      ir1_oo, ir2_oo, ir1_vv, ir2_vv, 
                      m_ham->m4get(i_oovv), m_int->m4get("iw_oooo"), m_int->m4get("iw_ooov"), m_int->m4get("iw2_ooov"), m_int->m4get("iw_ovov"), m_int->m4get("iw_ovvv"), m_int->m4get("iw2_ovvv"), m_int->m4get("iw_vvvv"));  

        a.set(r2, r2_new);
      }
    }
  }

  void character(std::vector<container_t> &v_rampl) const override {


    for (size_t ir = 0; ir < v_rampl.size(); ir++) {
      // normalise
      double norm = sqrt(v_rampl[ir].m2get(r1).dot(v_rampl[ir].m2get(r1)) + 0.25*v_rampl[ir].m4get(r2).dot(v_rampl[ir].m4get(r2)));
      v_rampl[ir].m2get(r1).scal(1/norm);
      v_rampl[ir].m4get(r2).scal(1/norm);
      double r12 = v_rampl[ir].m2get(r1).dot(v_rampl[ir].m2get(r1)) ;
      double r22 = 0.25*v_rampl[ir].m4get(r2).dot(v_rampl[ir].m4get(r2));

      std::ostringstream ss;

      // print excited state number and energy
      ss << "\n\nExcited state #" << ir+1 
         << "\n\nExcitation energy = " << std::setprecision(5) << std::fixed 
         << m_energy[ir] << " Ha = "
         << m_energy[ir]*inverse_electron_volt << " eV"
         << "\n\nr0 = "<< eom_ccsd_r0(m_energy[ir], v_rampl[ir].m2get(r1), v_rampl[ir].m4get(r2), m_int->m2get("if_ov"), m_ham->m4get(i_oovv)) 
         << "    ||r1||² = " << r12 << "    ||r2||² = " << r22
         << "\n\nAmplitude    Transition\n";
 
      std::vector<size_t> v_no, v_nv;
      gmb::get_tensor_dimensions(v_rampl[ir].m2get(r1), v_no, v_nv);

      // read r1  
      if ( r12 > r22 ) {
        constexpr size_t N = 2;
        libtensor::block_tensor_rd_ctrl<N, value_t> ctrl(v_rampl[ir].m2get(r1));

        libtensor::orbit_list<N, value_t> ol(ctrl.req_const_symmetry());
        for (libtensor::orbit_list<N, value_t>::iterator it = ol.begin(); it != ol.end(); it++) {
          libtensor::index<N> bidx;
          ol.get_index(it, bidx);
          libtensor::dense_tensor_rd_i<N, value_t> &blk = ctrl.req_const_block(bidx);
          libtensor::dense_tensor_rd_ctrl<N, value_t> tc(blk);
          const libtensor::dimensions<N> &tdims = blk.get_dims();
          const value_t *ptr = tc.req_const_dataptr();
          for (size_t offset = 0; offset < tdims.get_size(); offset++) {
            if (std::abs(ptr[offset]) >  0.1) {
              size_t i = 1+(offset/v_nv[bidx[1]]);
              size_t a = 1+offset-(offset/v_nv[bidx[1]])*v_nv[bidx[1]];
              ss << "\n" <<std::setw(8) << std::setprecision(5) << std::fixed <<  ptr[offset] << "     ";
              ss << "O" << i;
              for (size_t in = 0; in < N; in++) {
                ss << gmb::tospin(bidx[in]);
                if (bidx[in] > beta)
                  ss << bidx[in]-beta;
                if (in == 0)
                  ss << " -> V" << a;
              }
            }
          }
          tc.ret_const_dataptr(ptr);
          ctrl.ret_const_block(bidx);
        }
      } else {

      // read r2
      constexpr size_t N = 4;
      libtensor::block_tensor_rd_ctrl<N, value_t> ctrl(v_rampl[ir].m4get(r2));
      libtensor::orbit_list<N, value_t> ol(ctrl.req_const_symmetry());
      for (libtensor::orbit_list<N, value_t>::iterator it = ol.begin(); it != ol.end(); it++) {
        libtensor::index<N> bidx;
        ol.get_index(it, bidx);
        if (ctrl.req_is_zero_block(bidx)) 
          continue;
        libtensor::dense_tensor_rd_i<N, value_t> &blk = ctrl.req_const_block(bidx);
        libtensor::dense_tensor_rd_ctrl<N, value_t> tc(blk);
        const libtensor::dimensions<N> &tdims = blk.get_dims();
        const value_t *ptr = tc.req_const_dataptr();
        for (size_t offset = 0; offset < tdims.get_size(); offset++) {
          if (std::abs(ptr[offset]) >  0.1) {
            size_t i = offset / (v_no[bidx[1]]*v_nv[bidx[2]]*v_nv[bidx[3]]);
            size_t j = (offset - i*v_no[bidx[1]]*v_nv[bidx[2]]*v_nv[bidx[3]]) / (v_nv[bidx[2]]*v_nv[bidx[3]]);
            size_t a = (offset - j*v_nv[bidx[2]]*v_nv[bidx[3]] - i*v_no[bidx[1]]*v_nv[bidx[2]]*v_nv[bidx[3]]) / v_nv[bidx[3]];
            size_t b = offset - a*v_nv[bidx[3]] - j*v_nv[bidx[2]]*v_nv[bidx[3]] - i*v_no[bidx[1]]*v_nv[bidx[2]]*v_nv[bidx[3]];

            ss << "\n" << std::setw(8) << std::setprecision(5) << std::fixed <<  ptr[offset] << "     ";
            ss << "O" << 1+i << gmb::tospin(bidx[0]) << " -> V" << 1+a << gmb::tospin(bidx[2]);
            ss << "    O" << 1+j << gmb::tospin(bidx[1]) << " -> V" << 1+b << gmb::tospin(bidx[3]);
          }
        }
        tc.ret_const_dataptr(ptr);
        ctrl.ret_const_block(bidx);
      }
    }
    molpro::cout << ss.str() << "\n";
    }
  }

};

#endif // GMB_SRC_MOLPRO_GMB_PROBLEM_EOM_CCSD_H_