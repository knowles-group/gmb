#ifndef GMB_EOM_CCSD_IF_OO_H
#define GMB_EOM_CCSD_IF_OO_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief F_oo intermediate 
 * 
 * Needed to calculate CCSD T2s
 */
container<2,double> eom_ccsd_if_oo(
    container<2, double> &t1,     ///< CCSD T1
    container<4, double> &t2,     ///< CCSD T2
    // container<4, double> &tau,    ///< tau intermediate
    container<2, double> &f_oo,   ///< fock matrix - oo block
    container<2, double> &f_ov,   ///< fock matrix - ov block
    container<4, double> &i_ooov, ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv  ///> anti-symmetrized integral <ij||ab>
) {

    container<2, double> if_oo(f_oo.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    if_oo(i|j) = f_oo(i|j)
               + contract(e, t1(i|e), f_ov(j|e))    
               + contract(m|e, t1(m|e), i_ooov(j|m|i|e))      
               + contract(m|e|f, t1(i|e)*t1(m|f), i_oovv(j|m|e|f))
               + 0.5 * contract(m|e|f, t2(i|m|e|f), i_oovv(j|m|e|f))
               ;
    return if_oo;

}


#endif //GMB_EOM_CCSD_IF_OO_H
