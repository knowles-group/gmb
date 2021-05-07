#ifndef GMB_DIAG_OV_H
#define GMB_DIAG_OV_H

#include "../container.h"

/**
 * @brief Diagonal matrix containing energy differences.
 * 
 */

container<2,double> diag_ov(
    container<2, double> &f_oo,  ///< fock matrix - oo block
    container<2, double> &f_ov,  ///< fock matrix - ov block
    container<2, double> &f_vv  ///< fock matrix - vv block
) {
    container<2, double> d_ov(f_ov.get_space());

    libtensor::letter i,j,a,b;

    d_ov(i|a) = dirsum( diag(i, i|j, f_oo(i|j)), 
                         -diag(a, a|b, f_vv(a|b)));   

    return d_ov;
};

#endif // GMB_DIAG_OV_H
