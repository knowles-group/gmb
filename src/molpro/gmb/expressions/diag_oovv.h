#ifndef GMB_DIAG_OOVV_H
#define GMB_DIAG_OOVV_H

#include "../container.h"

/**
 * @brief Diagonal matrix containing energy differences.
 * 
 */

container<4,double> diag_oovv(
    container<2, double> &d_ov,    ///< diag - ov block
    container<4, double> &i_oovv   ///< anti-symmetrized integral
) {
    container<4, double> d_oovv(i_oovv.get_space());

    libtensor::letter i,j,a,b;

    d_oovv(i|j|a|b) = 0.25 * symm(a, b, symm(i, j, 
              dirsum(d_ov(i|a), d_ov(j|b))));

    return d_oovv;
};

#endif //GMB_DIAG_OOVV_H
