#ifndef PROBLEM_EOM_H_
#define PROBLEM_EOM_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "amplitudes.h"
#include "hamiltonian.h"
#include "expressions/ccsd/ccsd.h"


class problem_eom : public molpro::linalg::itsolv::Problem<supercontainer<>> {
protected:
  std::vector<double> m_energy; ///> energy
  mutable hamiltonian<> m_ham; ///> Hamiltonian
  mutable amplitudes<> m_tampl; ///> T amplitudes
public:
  using Problem::container_t;
  using Problem::value_t;
  problem_eom(const hamiltonian<> &ham, const amplitudes<> &tampl)
  : m_ham(ham), m_tampl(tampl) {}

  virtual ~problem_eom() {}

  void energy(container_t x) {
    // if (x.get_m2().find("t1") == x.get_m2().end())
    //   m_energy = 0.25 * x.m4get("t2").dot(m_ham.m4get("i_oovv"));
    // else
    //   m_energy = energy_ccsd(x.m2get("t1"), x.m4get("t2"), m_ham.m2get("f_ov"), m_ham.m4get("i_oovv"));
  }

  std::vector<double> get_energy() const {return m_energy;}  

  friend
  std::ostream& operator<<(std::ostream& s, const problem_eom& problem) ;
};

std::ostream& operator<<(std::ostream& s, const problem_eom& problem) {
  s<<"";
  return s;
}

#endif // PROBLEM_EOM_H_