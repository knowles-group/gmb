#ifndef GMB_SRC_MOLPRO_GMB_EOM_CCSD_IL_VV_H
#define GMB_SRC_MOLPRO_GMB_EOM_CCSD_IL_VV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief il_vv intermediate
 * 
 */
container<2,double> eom_ccsd_il_vv(
    container<2, double> &f_vv,     ///< fock operator - vv block
    container<4, double> &t2,       ///< CCSD T2 
    container<4, double> &l2        ///< EOM-CCSD L2 

) {

    container<2, double> il_vv(f_vv.get_space());

    libtensor::letter m,n,e,f,g;

    il_vv(e|f) = 0.5 * contract(m|n|g, l2(m|n|e|g), t2(m|n|f|g));
    return il_vv;

}

#endif // GMB_SRC_MOLPRO_GMB_EOM_CCSD_IL_VV_H
