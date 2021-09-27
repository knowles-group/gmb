#ifndef GMB_SRC_MOLPRO_GMB_D2_XXXX_H
#define GMB_SRC_MOLPRO_GMB_D2_XXXX_H

#include "../container.h"

/**
 * @brief Add dipole operator squared to 2-electron 
 * 
 */

void d2_xxxx(
    double fact,
    container<2, double> &d1_xx,
    container<2, double> &d2_xx,
    container<4, double> &int_xxxx
) {
    libtensor::letter p,q,r,s;
    int_xxxx(p|q|r|s) +=  fact*2.0*d1_xx(p|r)*b_xx(q|s) - a_xx(q|r)*a_xx(p|s);
}

#endif // GMB_SRC_MOLPRO_GMB_D2_XXXX_H
