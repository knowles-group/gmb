#ifndef ACTION_CCD_H
#define ACTION_CCD_H

#include <libtensor/libtensor.h>

/**
 * @brief Action to be applied to CCD T2 amplitudes.
 * 
 */

container<4,double> action_ccd (
    container<4, double> &t2,       ///< CCD T2 amplitude
    container<2, double> &f_oo,     ///< fock matrix - oo block
    container<2, double> &f_vv,     ///< fock matrix - vv block
    container<4, double> &i_oooo,   ///> anti-symmetrized integral <ij||kl>
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov,   ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_vvvv)   ///> anti-symmetrized integral <ab||cd>
{
    container<4, double> t2_new(t2.get_space()); 

    libtensor::letter i,j,m,n,a,b,e,f;

    t2_new(i|j|a|b) = i_oovv(i|j|a|b)
                    + asymm(a, b, contract(e, t2(i|j|a|e), f_vv(b|e)))
                    - asymm(i, j, contract(m, t2(i|m|a|b), f_oo(j|m)))
                    - asymm(a, b, asymm(i, j, 
                         contract(e|m, t2(i|m|a|e), i_ovov(m|b|j|e))))
                    + 0.5 * contract(m|n, t2(m|n|a|b), i_oooo(m|n|i|j))
                    + 0.5 * contract(e|f, t2(i|j|e|f), i_vvvv(a|b|e|f))
                    + 0.25 * contract(e|f, t2(i|j|e|f), 
                         contract(m|n, i_oovv(m|n|e|f), t2(m|n|a|b)))
                    + 0.5 * asymm(a, b, asymm(i, j, 
                         contract(m|e, t2(i|m|a|e), 
                             contract(n|f, i_oovv(m|n|e|f), t2(j|n|b|f)))))
                    - 0.5 * asymm(i, j, 
                         contract(n, t2(n|j|a|b), 
                             contract(m|e|f, i_oovv(m|n|e|f), t2(i|m|f|e))))
                    - 0.5 * asymm(a, b, 
                         contract(f, t2(i|j|f|b), 
                             contract(m|n|e, i_oovv(m|n|e|f), t2(n|m|a|e))));
 
    return t2_new;  
};


#endif // ACTION_CCD_H
