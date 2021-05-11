#ifndef GMB_EOM_CCSD_IR2_VV_H
#define GMB_EOM_CCSD_IR2_VV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief ir2_vv intermediate
 * 
 */
container<2,double> eom_ccsd_ir2_vv(
    container<2, double> &f_vv,     ///< fock operator - vv block
    container<4, double> &r2,       ///< EOM-CCSD R2
    container<4, double> &i_oovv    ///> anti-symmetrized integral <ij||ab>

) {

    container<2, double> ir2_vv(f_vv.get_space());

    libtensor::letter i,j,m,n,a,b,c,e,f;

    ir2_vv(a|e) = 0.5 * contract(m|n|f, r2(m|n|a|f), i_oovv(m|n|e|f));
    return ir2_vv;

}

#endif // GMB_EOM_CCSD_IR2_VV_H
