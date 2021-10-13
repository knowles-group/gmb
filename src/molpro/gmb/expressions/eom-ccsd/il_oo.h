#ifndef GMB_SRC_MOLPRO_GMB_EOM_CCSD_IL_OO_H
#define GMB_SRC_MOLPRO_GMB_EOM_CCSD_IL_OO_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief il_oo intermediate
 * 
 */
container<2,double> eom_ccsd_il_oo(
    container<2, double> &f_oo,     ///< fock operator - oo block
    container<4, double> &t2,       ///> CCSD T2
    container<4, double> &l2        ///< EOM-CCSD L2

) {

    container<2, double> il_oo(f_oo.get_space());

    libtensor::letter m,n,o,e,f;
    il_oo(m|n) = 0.5 * contract(e|f|o, l2(m|o|e|f), t2(n|o|e|f));

    return il_oo;

}

#endif // GMB_SRC_MOLPRO_GMB_EOM_CCSD_IL_OO_H
