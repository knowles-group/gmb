#ifndef GMB_PROBLEM_EOM_H_
#define GMB_PROBLEM_EOM_H_
#include <molpro/linalg/itsolv/IterativeSolver.h>
#include <vector>
#include "amplitudes.h"
#include "hamiltonian.h"
#include <molpro/iostream.h>

class problem_eom : public molpro::linalg::itsolv::Problem<amplitudes<>> {
protected:
  std::vector<double> m_energy;               ///> energy
  mutable hamiltonian<> m_ham;                ///> Hamiltonian
  mutable amplitudes<> m_tampl;               ///> T amplitudes
  size_t m_nroots;                            ///> energy
  mutable std::vector<amplitudes<>> m_vrampl; ///> T amplitudes
public:
  using Problem::container_t;
  using Problem::value_t;
  problem_eom(const hamiltonian<> &ham, const amplitudes<> &tampl, const size_t &nroots)
  : m_ham(ham), m_tampl(tampl), m_nroots(nroots) {
    m_vrampl.resize(nroots);
  }

  virtual ~problem_eom() {}

  void set_energy(std::vector<double> eigval) {m_energy = eigval;}
  std::vector<double> get_energy() const {return m_energy;}  

  /**
   * @brief Give description of transitions in terms of orbitals.
   * 
   */
  virtual void character() const {}

  };

#endif // GMB_PROBLEM_EOM_H_