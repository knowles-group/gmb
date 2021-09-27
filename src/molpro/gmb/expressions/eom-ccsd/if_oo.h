#ifndef GMB_SRC_MOLPRO_GMB_EOM_CCSD_IF_OO_H
#define GMB_SRC_MOLPRO_GMB_EOM_CCSD_IF_OO_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief F_oo intermediate 
 * 
 */
container<2,double> eom_ccsd_if_oo(
    container<2, double> &t1,     ///< CCSD T1
    container<4, double> &t2,     ///< CCSD T2
    container<2, double> &f_oo,   ///< fock matrix - oo block
    container<2, double> &f_ov,   ///< fock matrix - ov block
    container<4, double> &i_ooov, ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv  ///> anti-symmetrized integral <ij||ab>
) {

    container<2, double> if_oo(f_oo.get_space());

    libtensor::letter i,j,m,n,e,f;

    if_oo(i|j) = f_oo(i|j)
               + contract(e, t1(j|e), f_ov(i|e))    
               + contract(m|e, t1(m|e), i_ooov(i|m|j|e))      
               + contract(m|e|f, t1(j|e)*t1(m|f), i_oovv(i|m|e|f))
               + 0.5 * contract(m|e|f, t2(j|m|e|f), i_oovv(i|m|e|f));

    return if_oo;

}


#endif // GMB_SRC_MOLPRO_GMB_EOM_CCSD_IF_OO_H
