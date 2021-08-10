#ifndef GMB_PROBLEM_CCD_H_
#define GMB_PROBLEM_CCD_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "problem_gen.h"
#include "expressions/action_ccd.h"
#include "expressions/update.h"
#include "expressions/id_xxxx.h"


class problem_ccd : public problem_gen {
public:
  problem_ccd(const hamiltonian<> &ham)
  : problem_gen(ham) {}
  
  virtual ~problem_ccd() {}

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

  void action(const CVecRef<container_t> &parameters, const VecRef<container_t> &actions) const override {
    for (int k = 0; k < parameters.size(); k++) {
      auto &ccp = const_cast<container_t&> (parameters[k].get());     
      auto &a = actions[k].get();  
      auto t2_new = action_ccd(ccp.m4get(t2), m_ham.m2get(f_oo), m_ham.m2get(f_vv), 
                  m_ham.m4get(i_oooo), m_ham.m4get(i_oovv), m_ham.m4get(i_ovov), m_ham.m4get(i_vvvv));
      a.set(t2, t2_new);
    }
  }

  void energy(container_t x) {
      m_energy = 0.25 * x.m4get(t2).dot(m_ham.m4get(i_oovv));
  }

    void print(std::ostream& s) const {
    s << "CCD";
  }
};

#endif //GMB_PROBLEM_CCD_H_