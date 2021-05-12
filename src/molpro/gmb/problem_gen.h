#ifndef PROBLEM_GEN_H_
#define PROBLEM_GEN_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "hamiltonian.h"
#include "amplitudes.h"
#include "expressions/ccsd/energy.h"


class problem_gen : public molpro::linalg::itsolv::Problem<amplitudes<>> {
protected:
  mutable double m_energy = 0; ///> energy
  mutable hamiltonian<> m_ham; ///> Hamiltonian
public:
  using Problem::container_t;
  using Problem::value_t;
  problem_gen(const hamiltonian<> &ham)
  : m_ham(ham) {}

  virtual ~problem_gen() {}

  void energy(container_t x) {
    if (x.get_m2().find("t1") == x.get_m2().end())
      m_energy = 0.25 * x.m4get(t2).dot(m_ham.m4get(i_oovv));
    else
      m_energy = ccsd_energy(x.m2get(t1), x.m4get(t2), m_ham.m2get(f_ov), m_ham.m4get(i_oovv));
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