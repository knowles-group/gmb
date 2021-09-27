#ifndef GMB_SRC_MOLPRO_GMB_ADD_D2_H
#define GMB_SRC_MOLPRO_GMB_ADD_D2_H

#include "../container.h"

/**
 * @brief Add dipole operator squared to 2-electron 
 * 
 */

void add_d2(
    double fact,
    container<2, double> &d1_xx,
    container<2, double> &d2_xx,
    container<2, double> &d3_xx,
    container<2, double> &d4_xx,
    container<4, double> &int_xxxx
) {
    libtensor::letter p,q,r,s;
    int_xxxx(p|q|r|s) +=  fact*2.0*(d1_xx(p|r)*d2_xx(q|s) - d3_xx(q|r)*d4_xx(p|s));
}

#endif // GMB_SRC_MOLPRO_GMB_ADD_D2_H
