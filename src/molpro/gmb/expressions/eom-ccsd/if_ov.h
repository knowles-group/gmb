#ifndef GMB_EOM_CCSD_IF_OV_H
#define GMB_EOM_CCSD_IF_OV_H

#include <libtensor/libtensor.h>
#include "../../container.h"


/**
 * @brief F_ov intermediate
 * 
 */
container<2,double> eom_ccsd_if_ov(
    container<2, double> &t1,     ///< CCSD T1
    container<2, double> &f_ov,   ///< fock matrix - ov block
    container<4, double> &i_oovv  ///> anti-symmetrized integral <ij||ab>

) {

    container<2, double> if_ov(f_ov.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    if_ov(m|e) = f_ov(m|e)
                + contract(f|n, t1(n|f), i_oovv(m|n|e|f));

    return if_ov;

}


#endif //GMB_EOM_CCSD_IF_OV_H
