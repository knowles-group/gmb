#ifndef GMB_CCSD_IW_VVVV_H
#define GMB_CCSD_IW_VVVV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_vvvv intermediate
 * 
 */
container<4,double> ccsd_iw_vvvv(
    container<2, double> &t1,      ///< CCSD T1 amplitude
    container<4, double> &i_ovvv,  ///> anti-symmetrized integral <ia||bc>
    container<4, double> &i_vvvv   ///> anti-symmetrized integral <ab||cd>

) {

    container<4, double> iw_vvvv(i_vvvv.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    iw_vvvv(a|b|e|f) = i_vvvv(a|b|e|f)
                     - asymm(a, b, 
                       contract(m, t1(m|b), i_ovvv(m|a|f|e)));
                  
    return iw_vvvv;

}
#endif // GMB_CCSD_IW_VVVV_H
