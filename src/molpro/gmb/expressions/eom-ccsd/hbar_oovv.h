#ifndef GMB_HBAR_OOVV_H
#define GMB_HBAR_OOVV_H

#include "../../container.h"

/**
 * @brief Preconditioner for EOM - OOVV block
 * 
 */

container<4,double> hbar_oovv(
    container<2, double> &id_oo,   ///< identity matrix - oo block
    container<2, double> &id_vv,   ///< identity matrix - vv block
    container<4, double> &d_oovv,  ///< diag - oovv block
    double shift = 1e6             ///< shift
) {
    container<4, double> hb_oovv(d_oovv.get_space());

    libtensor::letter i,j,a,b;

    hb_oovv(i|j|a|b) = 
            0.25 * symm(i,j, symm(a,b,
            - d_oovv(i|j|a|b)
       	  + 0.5 * shift * symm(i, j, id_oo(i|j)) * symm(a, b, id_vv(a|b))));

    return hb_oovv;
};

#endif // GMB_HBAR_OOVV_H
