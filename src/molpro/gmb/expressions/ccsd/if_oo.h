#ifndef GMB_SRC_MOLPRO_GMB_CCSD_IF_OO_H
#define GMB_SRC_MOLPRO_GMB_CCSD_IF_OO_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief F_oo intermediate
 * 
 * Needed to calculate CCSD T1s
 * 
 */
container<2,double> ccsd_if_oo(
    container<2, double> &t1,     ///< CCSD T1
    container<4, double> &t2,     ///< CCSD T2 intermediate
    container<2, double> &f_oo,   ///< fock matrix - oo block
    container<4, double> &i_ooov, ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv  ///> anti-symmetrized integral <ij||ab>
) {

    container<2, double> if_oo(f_oo.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    if_oo(j|m) = f_oo(j|m)
               + contract(n|e, t1(n|e), i_ooov(j|n|m|e))      
               + 0.5 * contract(n|e|f, i_oovv(j|n|e|f), t2(m|n|e|f));
    return if_oo;
}


#endif // GMB_SRC_MOLPRO_GMB_CCSD_IF_OO_H
