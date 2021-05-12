#ifndef PROBLEM_CCD_H_
#define PROBLEM_CCD_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "problem_gen.h"
#include "action_ccd.h"
#include "update.h"


class problem_ccd : public problem_gen {
public:
  problem_ccd(const hamiltonian<> &ham)
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
    auto t2_new = action_ccd(ccv.m4get(t2), m_ham.m2get(f_oo), m_ham.m2get(f_vv), 
                  m_ham.m4get(i_oooo), m_ham.m4get(i_oovv), m_ham.m4get(i_ovov), m_ham.m4get(i_vvvv));
    a.set(t2, t2_new);
    return value;
  }

    void action(const CVecRef<container_t> &parameters, const VecRef<container_t> &actions) const override {}

  friend
  std::ostream& operator<<(std::ostream& s, const problem_ccd& problem) ;
};

std::ostream& operator<<(std::ostream& s, const problem_ccd& problem) {
  s<<"Problem: CCD ";
  return s;
}
#endif // PROBLEM_CCD_H_