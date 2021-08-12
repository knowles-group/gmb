#ifndef GMB_PROBLEM_EOM_CCSD_H_
#define GMB_PROBLEM_EOM_CCSD_H_
#include "problem_eom.h"
#include "expressions/eom-ccsd/eom-ccsd.h"
#include "expressions/update.h"

extern molpro::Profiler prof;
// extern std::unique_ptr<molpro::Profiler> pprof;

/**
 * @brief EOM-CCSD
 * 
 */

class problem_eom_ccsd : public problem_eom {
private:
  mutable supercontainer<> m_int;       ///> intermediates
public:
  problem_eom_ccsd(const hamiltonian<> &ham, const amplitudes<> &ampl)
  : problem_eom{ham, ampl} {
    init();
  }

  /**
   * @brief Initialise EOM calculation by computing EOM intermediates.
   * 
   */
  void init() {

    auto tau = ccsd_tau(m_tampl.m2get(t1), m_tampl.m4get(t2));
    auto if_oo = eom_ccsd_if_oo(m_tampl.m2get(t1), m_tampl.m4get(t2), m_ham.m2get(f_oo), m_ham.m2get(f_ov),
                   m_ham.m4get(i_ooov), m_ham.m4get(i_oovv));    
    auto if_ov = eom_ccsd_if_ov(m_tampl.m2get(t1), m_ham.m2get(f_ov), m_ham.m4get(i_oovv)); 
    auto if_vv = eom_ccsd_if_vv(m_tampl.m2get(t1), m_tampl.m4get(t2), m_ham.m2get(f_ov), m_ham.m2get(f_vv), 
                   m_ham.m4get(i_oovv), m_ham.m4get(i_ovvv));

    auto iw_oooo = eom_ccsd_iw_oooo(m_tampl.m2get(t1), tau, 
                    m_ham.m4get(i_oooo), m_ham.m4get(i_ooov), m_ham.m4get(i_oovv));    
    auto iw_ooov = eom_ccsd_iw_ooov(m_tampl.m2get(t1), m_tampl.m4get(t2), tau, if_ov,
                    iw_oooo, m_ham.m4get(i_ooov), m_ham.m4get(i_oovv), m_ham.m4get(i_ovov), m_ham.m4get(i_ovvv));
    auto iw_ovov = eom_ccsd_iw_ovov(m_tampl.m2get(t1), m_tampl.m4get(t2), 
                    m_ham.m4get(i_ooov), m_ham.m4get(i_oovv), m_ham.m4get(i_ovov), m_ham.m4get(i_ovvv));
    auto iw_vvvv = eom_ccsd_iw_vvvv(m_tampl.m2get(t1), tau, m_ham.m4get(i_oovv), m_ham.m4get(i_ovvv), m_ham.m4get(i_vvvv));
    auto iw_ovvv = eom_ccsd_iw_ovvv(m_tampl.m2get(t1), m_tampl.m4get(t2), tau, if_ov,  
                    m_ham.m4get(i_ooov), m_ham.m4get(i_oovv), m_ham.m4get(i_ovov), m_ham.m4get(i_ovvv), iw_vvvv);
    auto iw2_ooov = eom_ccsd_iw2_ooov(m_tampl.m2get(t1), m_ham.m4get(i_ooov), m_ham.m4get(i_oovv));
    auto iw2_ovvv = eom_ccsd_iw2_ovvv(m_tampl.m2get(t1), m_ham.m4get(i_oovv), m_ham.m4get(i_ovvv));
 
    m_int.set("tau", tau);
    m_int.set("if_oo", if_oo);
    m_int.set("if_ov", if_ov);
    m_int.set("if_vv", if_vv);
    m_int.set("iw_oooo", iw_oooo);
    m_int.set("iw_ooov", iw_ooov);
    m_int.set("iw_ovov", iw_ovov);
    m_int.set("iw_ovvv", iw_ovvv);
    m_int.set("iw_vvvv", iw_vvvv);
    m_int.set("iw2_ooov", iw2_ooov);
    m_int.set("iw2_ovvv", iw2_ovvv);   
  }

  bool diagonals(container_t &d) const override {
    auto d_ov = diag_ov(m_int.m2get("if_oo"), m_int.m2get("if_ov"), m_int.m2get("if_vv"));
    auto d_oovv = diag_oovv(d_ov, m_ham.m4get(i_oovv));
    auto id_oo = id_xx(m_ham.m2get(f_oo));
    auto id_vv = id_xx(m_ham.m2get(f_vv));
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
      auto id_oo = id_xx(m_ham.m2get(f_oo));
      auto id_vv = id_xx(m_ham.m2get(f_vv));
      auto a2 = precond_oovv(a.m4get(r2),d.m4get(r2), shift[k]);
      a.set(r1,a1);
      a.set(r2,a2);
    }
  }

  void action(const CVecRef<container_t> &parameters, const VecRef<container_t> &actions) const override {
    for (int k = 0; k < parameters.size(); k++) {
      auto &ccp = const_cast<container_t&> (parameters[k].get());     
      auto &a = actions[k].get();  
      // compute intermediates
      auto ir1_vv = eom_ccsd_ir1_vv(m_ham.m2get(f_vv),ccp.m2get(r1),m_int.m4get("iw2_ovvv"));          
      auto ir1_oo = eom_ccsd_ir1_oo(m_ham.m2get(f_oo),ccp.m2get(r1),m_int.m4get("iw2_ooov"));          
      auto ir2_oo = eom_ccsd_ir2_oo(m_ham.m2get(f_oo),ccp.m4get(r2),m_ham.m4get(i_oovv));          
      auto ir2_vv = eom_ccsd_ir2_vv(m_ham.m2get(f_vv),ccp.m4get(r2),m_ham.m4get(i_oovv));          
      // compute r1
      {
      auto r1_new = eom_ccsd_r1(ccp.m2get(r1), ccp.m4get(r2), m_int.m2get("if_oo"), m_int.m2get("if_ov"), m_int.m2get("if_vv"),  
                    m_int.m4get("iw_ovov"), m_int.m4get("iw2_ooov"), m_int.m4get("iw2_ovvv"));
      a.set(r1, r1_new);
      }
      // compute r2
      {
      auto r2_new = eom_ccsd_r2(ccp.m2get(r1), ccp.m4get(r2), m_tampl.m4get(t2), m_int.m2get("if_oo"), m_int.m2get("if_vv"),  
                    ir1_oo, ir2_oo, ir1_vv, ir2_vv, 
                    m_ham.m4get(i_oovv), m_int.m4get("iw_oooo"), m_int.m4get("iw_ooov"), m_int.m4get("iw2_ooov"), m_int.m4get("iw_ovov"), m_int.m4get("iw_ovvv"), m_int.m4get("iw2_ovvv"), m_int.m4get("iw_vvvv"));  
      
      a.set(r2, r2_new);
      }
    }
  }

    void character(std::vector<container_t> &v_rampl) const {
    constexpr size_t N =2;
    constexpr double inverse_electron_volt{27.211386245988};
    for (size_t ir1 = 0; ir1 < v_rampl.size(); ir1++) {
      // molpro::cout << "r[" << ir1 << "]:\n";
      // m_vr1[ir1]->print();
      molpro::cout << "\nExcited state #" << ir1+1 
                << "\nExcitation energy = " << std::setprecision(5) << std::fixed 
                << m_energy[ir1] << " Ha = "
                << m_energy[ir1]*inverse_electron_volt << " eV"
                << "\nocc -> vir     amplitude\n";

      libtensor::block_tensor_rd_i<2, value_t> &bt(v_rampl[ir1].m2get(r1));

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
        v_no.emplace_back(spl_o[i]-spl_o[i-1]);
    }
    v_no.emplace_back(no-std::accumulate(v_no.cbegin(),v_no.cend(),0));

    // virtual
    const libtensor::split_points &spl_v = bis.get_splits(1);
    for (size_t i = 0; i < spl_v.get_num_points(); i++) {
      if (i == 0)
        v_nv.push_back(spl_v[i]);
      else 
        v_nv.emplace_back(spl_v[i]-spl_v[i-1]);
    }
    v_nv.emplace_back(nv-std::accumulate(v_nv.cbegin(),v_nv.cend(),0));    

    std::vector<std::vector<size_t>> n_ne{v_no,v_nv};
      libtensor::block_tensor_rd_ctrl<N, value_t> ctrl(v_rampl[ir1].m2get(r1));

      libtensor::orbit_list<N, value_t> ol(ctrl.req_const_symmetry());
      for (libtensor::orbit_list<N, value_t>::iterator it = ol.begin(); it != ol.end(); it++) {
        libtensor::index<N> bidx;
        ol.get_index(it, bidx);
        const libtensor::dimensions<N> &bd = bt.get_bis().get_block_dims(bidx);
        libtensor::dense_tensor_rd_i<N, value_t> &blk = ctrl.req_const_block(bidx);
        libtensor::dense_tensor_rd_ctrl<N, value_t> tc(blk);
        const libtensor::dimensions<N> &tdims = blk.get_dims();
        const value_t *ptr = tc.req_const_dataptr();
        for (size_t offset = 0; offset < tdims.get_size(); offset++) {
          if (std::abs(ptr[offset]) >  0.001) {
            size_t i = 1+(offset/v_nv[bidx[1]]);
            size_t a = 1+offset-(offset/v_nv[bidx[1]])*v_nv[bidx[1]];
            molpro::cout << "o" << i;
            for (size_t in = 0; in < N; in++) {
              switch (bidx[in]) {
              case alpha: molpro::cout << "a";
                break;
              case beta: molpro::cout << "b";
                break;
              default: molpro::cout << "p" << bidx[in]-beta;
                break;
              }
              if (in == 0)
                molpro::cout << " -> v" << a;
            }
            molpro::cout << "     " << std::setprecision(5) << std::fixed <<  ptr[offset] << "\n";
          }
        }
        tc.ret_const_dataptr(ptr);
        ctrl.ret_const_block(bidx);
      }
    }
  };


};

#endif //GMB_PROBLEM_EOM_CCSD_H_