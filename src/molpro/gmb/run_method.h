#ifndef RUN_METHOD_H
#define RUN_METHOD_H

#include "hamiltonian.h"
#include "amplitudes.h"
#include "problem_gen.h"
#include "problem_eom.h"

/**
 * @brief Run Ground-State (GS) calculation
 * 
 */
void run_gs(hamiltonian<> &ham, const std::string &method, std::unique_ptr<problem_gen> &problem, std::unique_ptr<amplitudes<>> &ptampl);

/**
 * @brief Run Excited-State (ES) calculation
 * 
 */
void run_es(const hamiltonian<> &ham, const std::string &method, std::unique_ptr<problem_eom> &problem, const std::unique_ptr<amplitudes<>> &ptampl, const size_t &nroots);

#endif // RUN_METHOD_H