#ifndef GMB_SRC_MOLPRO_GMB_PROBLEM_EOM_CCSD_LEFT_H_
#define GMB_SRC_MOLPRO_GMB_PROBLEM_EOM_CCSD_LEFT_H_
#include "constants.h"
#include "expressions/eom-ccsd/eom-ccsd-left.h"
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

class problem_eom_ccsd_left : public problem_eom {
private:
  const std::shared_ptr<supercontainer<>> m_int;       ///> intermediates
public:
  problem_eom_ccsd_left(const std::shared_ptr<hamiltonian<>> &ham, 
                   const std::shared_ptr<amplitudes<>> &ptampl,
                   const std::shared_ptr<supercontainer<>> &intermediates)
  : problem_eom{ham, ptampl}, m_int(intermediates) {
  }

  void action(const CVecRef<container_t> &parameters, const VecRef<container_t> &actions) const override {
  
    for (int k = 0; k < parameters.size(); k++) {

      auto &ccp = const_cast<container_t&> (parameters[k].get());     
      auto &a = actions[k].get();  

      singlet_projector(ccp);

      // compute left intermediates
      auto il_vv = eom_ccsd_il_vv(m_ham->m2get(f_vv), m_tampl->m4get(t2), ccp.m4get(r2));          
      auto il_oo = eom_ccsd_il_oo(m_ham->m2get(f_oo), m_tampl->m4get(t2), ccp.m4get(r2));          

      // compute l1
      {
        auto l1_new = eom_ccsd_l1(ccp.m2get(r1), 
                                  ccp.m4get(r2), 
                                  m_tampl->m4get(t2), 
                                  il_oo, il_vv, 
                                  m_int->m2get("if_oo"), 
                                  m_int->m2get("if_ov"), 
                                  m_int->m2get("if_vv"),  
                                  m_int->m4get("iw_ovov"), 
                                  m_int->m4get("iw_ooov"), 
                                  m_int->m4get("iw2_ooov"), 
                                  m_int->m4get("iw_ovvv"), 
                                  m_int->m4get("iw2_ovvv"));
        a.set(r1, l1_new);
      }
      // compute l2
      {
        auto l2_new = eom_ccsd_l2(ccp.m2get(r1), 
                                  ccp.m4get(r2), 
                                  m_tampl->m4get(t2), 
                                  m_int->m2get("if_oo"), 
                                  m_int->m2get("if_ov"), 
                                  m_int->m2get("if_vv"),  
                                  il_oo, il_vv, 
                                  m_ham->m4get(i_oovv), 
                                  m_int->m4get("iw_oooo"), 
                                  m_int->m4get("iw_ooov"), 
                                  m_int->m4get("iw2_ooov"), 
                                  m_int->m4get("iw_ovov"), 
                                  m_int->m4get("iw_ovvv"), 
                                  m_int->m4get("iw2_ovvv"), 
                                  m_int->m4get("iw_vvvv"));  
        a.set(r2, l2_new);
      }
    }
  }


};

#endif // GMB_SRC_MOLPRO_GMB_PROBLEM_EOM_CCSD_LEFT_H_