#ifndef PROBLEM_CCSD_H_
#define PROBLEM_CCSD_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "problem_gen.h"
#include "action_ccsd.h"
#include "update.h"


class problem_ccsd : public problem_gen {
public:
  problem_ccsd(const supercontainer<> &ham)
  : problem_gen(ham) {}

  void precondition(const VecRef<container_t> &g,
                    const std::vector<value_t> &shift) const override {
    for (int k = 0; k < g.size(); k++) {
      auto &a = g[k].get();     
      auto t1_new = update_t1((m_ham.m2get("f_oo")), (m_ham.m2get("f_vv")), a.m2get("t1"));
      auto t2_new = update_t2(m_ham.m2get("f_oo"), m_ham.m2get("f_vv"), a.m4get("t2"));
      a.set("t1", t1_new);
      a.set("t2", t2_new);
    }
  }


  value_t residual(const container_t &v, container_t &a) const override {
    value_t value = 0;
    auto &ccv = const_cast<container_t&> (v);    
    auto t1_new = action_ccsd_t1(ccv.m2get("t1"), ccv.m4get("t2"), m_ham.m2get("f_oo"), m_ham.m2get("f_ov"), m_ham.m2get("f_vv"), 
                   m_ham.m4get("i_oooo"), m_ham.m4get("i_ooov"), m_ham.m4get("i_oovv"), m_ham.m4get("i_ovov"), 
                   m_ham.m4get("i_ovvv"), m_ham.m4get("i_vvvv"));
    auto t2_new = action_ccsd_t2(ccv.m2get("t1"), ccv.m4get("t2"), m_ham.m2get("f_oo"), m_ham.m2get("f_ov"), m_ham.m2get("f_vv"), 
                    m_ham.m4get("i_oooo"), m_ham.m4get("i_ooov"), m_ham.m4get("i_oovv"), m_ham.m4get("i_ovov"), 
                    m_ham.m4get("i_ovvv"), m_ham.m4get("i_vvvv"));
    a.set("t1", t1_new);
    a.set("t2", t2_new);  
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