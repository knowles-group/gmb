#ifndef GMB_HBAR_OV_H
#define GMB_HBAR_OV_H

#include "../../container.h"

/**
 * @brief \bar{H} EOM - OOVV block
 * 
 */

container<2,double> hbar_ov(
    container<2, double> &d_ov  ///< diag - ov block
) {
    container<2, double> hb_ov(d_ov.get_space());

    libtensor::letter i,j,a,b;

    hb_ov(i|a) = - d_ov(i|a);

    return hb_ov;
};

#endif // GMB_HBAR_OV_H
