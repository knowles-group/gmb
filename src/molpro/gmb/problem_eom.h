#ifndef GMB_SRC_MOLPRO_GMB_PROBLEM_EOM_H_
#define GMB_SRC_MOLPRO_GMB_PROBLEM_EOM_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <molpro/iostream.h>
#include "amplitudes.h"
#include "hamiltonian.h"

class problem_eom : public molpro::linalg::itsolv::Problem<amplitudes<>> {

protected:
  double m_e0{0.0};
  std::vector<value_t> m_energies; ///> energy
  mutable hamiltonian<> m_ham;   ///> Hamiltonian
  mutable amplitudes<> m_tampl;  ///> T amplitudes

public:
  using Problem::container_t;
  using Problem::value_t;

  problem_eom(const hamiltonian<> &ham, const amplitudes<> &tampl)
  : m_ham{ham}, m_tampl{tampl} {
  }

  virtual ~problem_eom() = default;

  void set_e0(double e0) {m_e0 = e0;}

  void set_energy(std::vector<value_t> eigval) {m_energies = eigval;}
  std::vector<value_t> get_energy() const {return m_energies;}  

  /**
   * @brief Give description of transitions in terms of orbitals
   * 
   */
  virtual void character(std::vector<container_t> &v_rampl) const {}

  };

#endif // GMB_PROBLEM_EOM_H_