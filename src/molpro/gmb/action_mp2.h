#ifndef ACTION_MP2_H
#define ACTION_MP2_H

#include <libtensor/libtensor.h>

/**
 * @brief Action to be applied to MP2 T2 amplitudes.
 * 
 */

container<4,double> action_mp2 (
    container<4, double> &t2,     ///< MP2 T2 amplitude
    container<2, double> &f_oo,  ///< fock matrix - oo block
    container<2, double> &f_vv   ///< fock matrix - vv block
) {
    container<4, double> t2_new(t2.get_space());

    libtensor::letter i,j,m,a,b,e;
    t2_new(i|j|a|b) = asymm(a, b, contract(e, t2(i|j|a|e), f_vv(b|e)))
                    - asymm(i, j, contract(m, t2(i|m|a|b), f_oo(m|j)));
    return t2_new;  
};

#endif // ACTION_MP2_H
