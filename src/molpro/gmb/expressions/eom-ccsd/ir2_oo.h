#ifndef GMB_EOM_CCSD_IR2_OO_H
#define GMB_EOM_CCSD_IR2_OO_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief ir2_oo intermediate
 * 
 */
container<2,double> eom_ccsd_ir2_oo(
    container<2, double> &f_oo,     ///< fock operator - vv block
    container<4, double> &r2,       ///< EOM-CCSD R2
    container<4, double> &i_oovv    ///> anti-symmetrized integral <ij||ab>

) {

    container<2, double> ir2_oo(f_oo.get_space());

    libtensor::letter i,j,m,n,a,b,c,e,f;
    ir2_oo(j|m) = 0.5 * contract(n|e|f, r2(j|n|e|f), i_oovv(m|n|e|f));

    return ir2_oo;

}

#endif //GMB_EOM_CCSD_IR2_OO_H
