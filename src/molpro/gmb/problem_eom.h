#ifndef GMB_PROBLEM_EOM_H_
#define GMB_PROBLEM_EOM_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/iostream.h>
#include "amplitudes.h"
#include "hamiltonian.h"

class problem_eom : public molpro::linalg::itsolv::Problem<amplitudes<>> {

protected:
  std::vector<value_t> m_energy; ///> energy
  mutable hamiltonian<> m_ham;   ///> Hamiltonian
  mutable amplitudes<> m_tampl;  ///> T amplitudes

public:
  using Problem::container_t;
  using Problem::value_t;

  problem_eom(const hamiltonian<> &ham, const amplitudes<> &tampl)
  : m_ham{ham}, m_tampl{tampl} {}

  virtual ~problem_eom() {}

  void set_energy(std::vector<value_t> eigval) {m_energy = eigval;}
  std::vector<value_t> get_energy() const {return m_energy;}  

  /**
   * @brief Give description of transitions in terms of orbitals
   * 
   */
  virtual void character(std::vector<container_t> &v_rampl) const {}

  virtual void check_eigenvalue(const container_t &rampl) const {}

  };

#endif // GMB_PROBLEM_EOM_H_