#ifndef GMB_PROBLEM_MP2_H_
#define GMB_PROBLEM_MP2_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "problem_gen.h"
#include "action_mp2.h"
#include "update.h"
#include "expressions/id_xxxx.h"


class problem_mp2 : public problem_gen {
public:
  problem_mp2(const hamiltonian<> &ham)
  : problem_gen(ham) {}

  void precondition(const VecRef<container_t> &g,
                    const std::vector<value_t> &shift) const override {
    for (int k = 0; k < g.size(); k++) {
      auto &a = g[k].get();     
      auto t2_new = update_t2(m_ham.m2get(f_oo), m_ham.m2get(f_vv), a.m4get(t2));
      a.set(t2, t2_new);
    }
  }


  value_t residual(const container_t &v, container_t &a) const override {
    value_t value = 0;
    auto &ccv = const_cast<container_t&> (v);    
    auto t2_new = action_mp2(ccv.m4get(t2), m_ham.m2get(f_oo), m_ham.m2get(f_vv));
    a.set(t2, t2_new);
    return value;
  }

  void action(const CVecRef<container_t> &parameters, const VecRef<container_t> &actions) const override {
    for (int k = 0; k < parameters.size(); k++) {
      auto &ccp = const_cast<container_t&> (parameters[k].get());     
      auto &a = actions[k].get();  
      auto t2_new = action_mp2(ccp.m4get(t2), m_ham.m2get(f_oo), m_ham.m2get(f_vv), 
                  m_ham.m4get(i_oooo), m_ham.m4get(i_oovv), m_ham.m4get(i_ovov), m_ham.m4get(i_vvvv));
      a.set(t2, t2_new);
    }
  }

  friend
  std::ostream& operator<<(std::ostream& s, const problem_mp2& problem) ;
};

std::ostream& operator<<(std::ostream& s, const problem_mp2& problem) {
  s<<"Problem: MP2 ";
  return s;
}
#endif //GMB_PROBLEM_MP2_H_