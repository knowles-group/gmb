#ifndef PROBLEM_CCSD_H_
#define PROBLEM_CCSD_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "problem_gen.h"
#include "expressions/ccsd/ccsd.h"
#include "update.h"


class problem_ccsd : public problem_gen {
public:
  problem_ccsd(const hamiltonian<> &ham)
  : problem_gen(ham) {}

  void precondition(const VecRef<container_t> &g,
                    const std::vector<value_t> &shift) const override {
    for (int k = 0; k < g.size(); k++) {
      auto &a = g[k].get();     
      auto t1_new = update_t1((m_ham.m2get(f_oo)), (m_ham.m2get(f_vv)), a.m2get(t1));
      auto t2_new = update_t2(m_ham.m2get(f_oo), m_ham.m2get(f_vv), a.m4get(t2));
      a.set(t1, t1_new);
      a.set(t2, t2_new);
    }
  }


  value_t residual(const container_t &v, container_t &a) const override {
    value_t value = 0;
    auto &ccv = const_cast<container_t&> (v);        

    auto if_oo = ccsd_if_oo(ccv.m2get(t1), ccv.m4get(t2), m_ham.m2get(f_oo),
                   m_ham.m4get(i_ooov), m_ham.m4get(i_oovv));    
    auto if_ov = ccsd_if_ov(ccv.m2get(t1), m_ham.m2get(f_ov), m_ham.m4get(i_oovv)); 
    auto if_vv = ccsd_if_vv(ccv.m2get(t1), ccv.m4get(t2), m_ham.m2get(f_ov), m_ham.m2get(f_vv), 
                   m_ham.m4get(i_oovv), m_ham.m4get(i_ovvv));

    auto if2_oo = ccsd_if2_oo(ccv.m2get(t1), ccv.m4get(t2), m_ham.m2get(f_oo), m_ham.m2get(f_ov),
                    m_ham.m4get(i_ooov), m_ham.m4get(i_oovv));    
    auto tau = ccsd_tau(ccv.m2get(t1), ccv.m4get(t2));
    auto iw_oooo = ccsd_iw_oooo(ccv.m2get(t1), tau,   
                    m_ham.m4get(i_oooo), m_ham.m4get(i_ooov), m_ham.m4get(i_oovv));
    auto iw_vvvv = ccsd_iw_vvvv(ccv.m2get(t1), 
                    m_ham.m4get(i_ovvv), m_ham.m4get(i_vvvv));
    auto iw_ovov = ccsd_iw_ovov(ccv.m2get(t1), ccv.m4get(t2), 
                    m_ham.m4get(i_ooov), m_ham.m4get(i_oovv), m_ham.m4get(i_ovov), m_ham.m4get(i_ovvv));
                         
    auto t1_new = ccsd_t1(ccv.m2get(t1), ccv.m4get(t2), m_ham.m2get(f_ov), if_oo, if_ov, if_vv, 
                    m_ham.m4get(i_ooov),m_ham.m4get(i_ovov), m_ham.m4get(i_ovvv));
    auto t2_new = ccsd_t2(ccv.m2get(t1), ccv.m4get(t2), tau, if2_oo, if_vv, 
                    m_ham.m4get(i_ooov), m_ham.m4get(i_oovv), m_ham.m4get(i_ovov), m_ham.m4get(i_ovvv), 
                    iw_oooo, iw_ovov, iw_vvvv);
           
    a.set(t1, t1_new);
    a.set(t2, t2_new);  
    return value;
  }

    void action(const CVecRef<container_t> &parameters, const VecRef<container_t> &actions) const override {}

  friend
  std::ostream& operator<<(std::ostream& s, const problem_ccsd& problem) ;
};

std::ostream& operator<<(std::ostream& s, const problem_ccsd& problem) {
  s<<"Problem: CCSD ";
  return s;
}
#endif // PROBLEM_CCSD_H_
