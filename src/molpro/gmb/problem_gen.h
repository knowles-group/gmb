#ifndef PROBLEM_GEN_H_
#define PROBLEM_GEN_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "supercontainer.h"
#include "action_ccsd.h"


class problem_gen : public molpro::linalg::itsolv::Problem<supercontainer<>> {
protected:
  mutable double m_energy = 0; ///> energy
  mutable supercontainer<> m_ham; ///> Hamiltonian
public:
  using Problem::container_t;
  using Problem::value_t;
  problem_gen(const supercontainer<> &ham)
  : m_ham(ham) {}

  void energy(container_t x) {
    if (x.get_m2().find("t1") == x.get_m2().end())
      m_energy = 0.25 * x.m4get("t2").dot(m_ham.m4get("i_oovv"));
    else
      m_energy = energy_ccsd(x.m2get("t1"), x.m4get("t2"), m_ham.m2get("f_ov"), m_ham.m4get("i_oovv"));
  }

  double get_energy() const {return m_energy;}  

  friend
  std::ostream& operator<<(std::ostream& s, const problem_gen& problem) ;
};

std::ostream& operator<<(std::ostream& s, const problem_gen& problem) {
  s<<"";
  return s;
}

#endif // PROBLEM_GEN_H_