#ifndef GMB_SRC_MOLPRO_GMB_RUN_METHOD_H
#define GMB_SRC_MOLPRO_GMB_RUN_METHOD_H

#include "hamiltonian.h"
#include "amplitudes.h"
#include "problem_gen.h"
#include "problem_eom.h"

/**
 * @brief Run Ground-State (GS) calculation
 * 
 */
void run_gs(hamiltonian<> &ham, 
            const std::string &method, 
            std::unique_ptr<problem_gen> &problem, 
            std::shared_ptr<amplitudes<>> &ptampl);

/**
 * @brief Run equation-of-motion (EOM) coupled-cluster (CC) calculation
 * 
 */
void run_eom(const std::shared_ptr<hamiltonian<>> &pham, 
             const std::string &method, 
             std::unique_ptr<problem_eom> &problem, 
             const std::shared_ptr<amplitudes<>> &ptampl, 
             const size_t &nroots, 
             const double& es_conv);

/**
 * @brief Compute EOM intermediates.
 * 
 */
void eom_intermediates(const std::shared_ptr<hamiltonian<>> &pham, 
                       const std::shared_ptr<amplitudes<>> &tampl, 
                       std::shared_ptr<supercontainer<>>  &intermediates);

#endif // GMB_SRC_MOLPRO_GMB_RUN_METHOD_H