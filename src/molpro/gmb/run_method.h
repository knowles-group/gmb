#ifndef RUN_METHOD_H
#define RUN_METHOD_H

#include "hamiltonian.h"
#include "amplitudes.h"
#include "problem_gen.h"
#include "problem_eom.h"

void run_cc(hamiltonian<> &ham, const std::string &method, std::unique_ptr<problem_gen> &problem, std::unique_ptr<amplitudes<>> &ptampl);

void run_eom(const hamiltonian<> &ham, const std::string &method, std::unique_ptr<problem_eom> &problem, const std::unique_ptr<amplitudes<>> &ptampl, const size_t &nroots);

#endif // RUN_METHOD_H