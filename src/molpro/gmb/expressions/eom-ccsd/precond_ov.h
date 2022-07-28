#ifndef GMB_SRC_MOLPRO_GMB_PRECOND_OV_H
#define GMB_SRC_MOLPRO_GMB_PRECOND_OV_H

#include "../../container.h"

/**
 * @brief Preconditioner for EOM - OOVV block
 * 
 */

container<2,double> precond_ov(
    container<2, double> &res_ov,  ///< residual - ov block
    container<2, double> &d_ov,    ///< diag - ov block
    double shift                   ///< shift
) {
    container<2, double> pc_ov(d_ov.get_space());

    libtensor::letter i,j,a,b;

    pc_ov(i|a) = - div(res_ov(i|a), 
                     (d_ov(i|a) - set((double)shift+1e-12, d_ov(i|a))));

    return pc_ov;
};

#endif // GMB_SRC_MOLPRO_GMB_PRECOND_OV_H
