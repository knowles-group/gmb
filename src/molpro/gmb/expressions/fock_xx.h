#ifndef FOCK_XX_H
#define FOCK_XX_H

#include "../container.h"

/**
 * @brief Compute Fock matrix - xx block
 * 
 */

container<2,double> fock_xx (
    container<2, double> &d_oo,    ///> identity matrix
    container<2, double> &h1_xx,   ///> one-electron hamiltonian
    container<4, double> &i_oxox   ///> anti-symmetrized integral <ip||jq>
) {
    container<2,double> f_xx(h1_xx.get_space());
    libtensor::letter i,j,p,q;
    f_xx(p|q) = h1_xx(p|q)
              + contract(i|j, d_oo(i|j), i_oxox(i|p|j|q));
    return f_xx;
}

#endif // FOCK_XX_H
