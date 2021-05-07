#ifndef GMB_PRECOND_OOVV_H
#define GMB_PRECOND_OOVV_H

#include "../../container.h"

/**
 * @brief Preconditioner for EOM - OOVV block
 * 
 */

container<4,double> precond_oovv(
    container<4, double> &res_oovv,  ///< diag - oovv block
    container<4, double> &d_oovv,    ///< diag - oovv block
    double shift                     ///< shift
) {
    container<4, double> pc_oovv(d_oovv.get_space());

    libtensor::letter i,j,a,b;

    pc_oovv(i|j|a|b) =  div(res_oovv(i|j|a|b), 
                     (d_oovv(i|j|a|b) - set((double)shift, d_oovv(i|j|a|b))));
    return pc_oovv;
};

#endif // GMB_PRECOND_OOVV_H
