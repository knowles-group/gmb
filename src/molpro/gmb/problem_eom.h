#ifndef GMB_PROBLEM_EOM_H_
#define GMB_PROBLEM_EOM_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "amplitudes.h"
#include "hamiltonian.h"
#include "expressions/ccsd/ccsd.h"


class problem_eom : public molpro::linalg::itsolv::Problem<amplitudes<>> {
protected:
  std::vector<double> m_energy;    ///> energy
  mutable hamiltonian<> m_ham;     ///> Hamiltonian
  mutable amplitudes<> m_tampl;    ///> T amplitudes
public:
  using Problem::container_t;
  using Problem::value_t;
  problem_eom(const hamiltonian<> &ham, const amplitudes<> &tampl)
  : m_ham(ham), m_tampl(tampl) {}

  virtual ~problem_eom() {}

  void set_energy(std::vector<double> eigval) {m_energy = eigval;}

  std::vector<double> get_energy() const {return m_energy;}  

  friend
  std::ostream& operator<<(std::ostream& s, const problem_eom& problem) ;

  virtual 
  void print(std::ostream& s) const {}

};

std::ostream& operator<<(std::ostream& s, const problem_eom& problem) {
  problem.print(s);
  return s;
}

#endif // GMB_PROBLEM_EOM_H_