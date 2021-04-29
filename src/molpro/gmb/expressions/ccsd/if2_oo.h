#ifndef GMB_CCSD_IF2_OO_H
#define GMB_CCSD_IF2_OO_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief F_oo intermediate - for t2
 * 
 * Needed to calculate CCSD T2 amplitudes
 */
container<2,double> ccsd_if2_oo(
    container<2, double> &t1,     ///< CCSD T1 amplitude
    container<4, double> &t2,     ///< CCSD T2 intermediate
    container<2, double> &f_oo,   ///< fock matrix - oo block
    container<2, double> &f_ov,   ///< fock matrix - ov block
    container<4, double> &i_ooov, ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv  ///> anti-symmetrized integral <ij||ab>
) {

    container<2, double> if_oo(f_oo.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    if_oo(j|m) = f_oo(j|m)
               + contract(n|e, t1(n|e), i_ooov(m|n|j|e))      
               + contract(e, t1(j|e), f_ov(m|e))    
               + 0.5 * contract(n|e|f, i_oovv(n|m|e|f), t2(j|n|f|e))
               + contract(f|n, t1(n|f), 
                    contract(e, t1(j|e), i_oovv(m|n|e|f)));
    return if_oo;

}


#endif // GMB_CCSD_IF2_OO_H
