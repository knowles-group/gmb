#ifndef GMB_CCSD_T1_H
#define GMB_CCSD_T1_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief Action to be applied to CCSD T1 amplitudes.
 * 
 */
container<2,double> ccsd_t1(
    container<2, double> &t1,      ///< CCSD T1 amplitude
    container<4, double> &t2,      ///< CCSD T2 amplitude
    container<2, double> &f_ov,    ///< fock matrix - ov block
    container<2, double> &if_oo,   ///< intermediate - oo block
    container<2, double> &if_ov,   ///< intermediate - ov block
    container<2, double> &if_vv,   ///< intermediate - vv block
    container<4, double> &i_ooov,  ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_ovov,  ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_ovvv   ///> anti-symmetrized integral <ia||bc>

) {

    container<2, double> t1_new(t1.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    t1_new(i|a) = f_ov(i|a)
                + contract(e, t1(i|e), if_vv(a|e))
                - contract(m, t1(m|a), if_oo(m|i))
                + contract(e|m, t2(i|m|a|e), if_ov(m|e))
                - contract(e|m, t1(m|e), i_ovov(m|a|i|e))  
                - 0.5 * contract(e|f|m, t2(i|m|e|f), i_ovvv(m|a|e|f))
                - 0.5 * contract(e|m|n, t2(m|n|a|e), i_ooov(m|n|i|e));

    return t1_new;

}

#endif // GMB_CCSD_T1_H
