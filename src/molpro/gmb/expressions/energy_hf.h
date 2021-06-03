#ifndef ENERGY_HF_H
#define ENERGY_HF_H

#include "../container.h"

/**
 * @brief Compute Hartree-Fock (HF) energy
 * 
 */

double energy_hf(
    container<2, double> &f_oo, ///> Fock matrix - oo block
    container<4, double> &i_oooo ///> anti-symmetrized integral <ij||kl>
) {

    double energy;

    libtensor::letter i, j, k, l;
    energy = trace(i,j,f_oo(i|j))
                  - 0.5*trace (i|j, k|l, i_oooo(i|j|k|l))
                  ;
    return energy;
};

#endif // ENERGY_HF_H
