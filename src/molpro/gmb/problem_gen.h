#ifndef GMB_PROBLEM_GEN_H_
#define GMB_PROBLEM_GEN_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "hamiltonian.h"
#include "amplitudes.h"

class problem_gen : public molpro::linalg::itsolv::Problem<amplitudes<>> {
protected:
  double m_energy{0.0}; ///> energy
  mutable hamiltonian<> m_ham; ///> Hamiltonian
public:
  using Problem::container_t;
  using Problem::value_t;
  problem_gen(const hamiltonian<> &ham)
  : m_ham{ham} {}

  virtual ~problem_gen() = default;


  virtual void energy(container_t x) {};

  double get_energy() const {
    return m_energy;
  }  

  friend
  std::ostream& operator<<(std::ostream& s, const problem_gen& problem);

  virtual 
  void print(std::ostream& s) const {}
};

// std::ostream& operator<<(std::ostream& s, const problem_gen& problem) {
//   problem.print(s);
//   return s;
// }

#endif //GMB_PROBLEM_GEN_H_