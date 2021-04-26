#ifndef ACTION_CCSD_H
#define ACTION_CCSD_H

#include <libtensor/libtensor.h>
#include "container.h"


double energy_ccsd (
    libtensor::any_tensor<2, double> &t1,     ///< CCSD T1 amplitude
    libtensor::any_tensor<4, double> &t2,     ///< CCSD T2 amplitude
    libtensor::any_tensor<2, double> &f_ov,   ///< fock matrix - ov block
    libtensor::any_tensor<4, double> &i_oovv  ///> anti-symmetrized integral <ij||ab>
) {
    
    libtensor::letter i,j,m,n,a,b,e,f;

    double energy = dot_product(t1(m|e), f_ov(m|e))
               + 0.5 * dot_product(t1(m|e), 
                    contract(n|f, t1(n|f), i_oovv(m|n|e|f)))
               + 0.25 * dot_product(t2(m|n|e|f), i_oovv(m|n|e|f));
    
   return energy;
}


/**
 * @brief Action to be applied to CCSD T1 amplitudes.
 * 
 */
container<2,double> action_ccsd_t1(
    container<2, double> &t1,     ///< CCSD T1 amplitude
    container<4, double> &t2,     ///< CCSD T2 amplitude
    container<2, double> &f_oo,   ///< fock matrix - oo block
    container<2, double> &f_ov,   ///< fock matrix - ov block
    container<2, double> &f_vv,   ///< fock matrix - vv block
    container<4, double> &i_oooo, ///> anti-symmetrized integral <ij||kl>
    container<4, double> &i_ooov, ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv, ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov, ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_ovvv, ///> anti-symmetrized integral <ia||bc>
    container<4, double> &i_vvvv  ///> anti-symmetrized integral <ab||cd>

) {

    container<2, double> t1_new(t1.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    t1_new(i|a) = f_ov(i|a)
                   - contract(m, t1(m|a), f_oo(i|m))
                   + contract(e, t1(i|e), f_vv(a|e))
                   - contract(e|m, t1(m|e), i_ovov(m|a|i|e))
                   + contract(e|m, t2(i|m|a|e), f_ov(m|e))
                   - 0.5 * contract(e|m|n, t2(m|n|a|e), i_ooov(m|n|i|e))
                   + 0.5 * contract(e|f|m, t2(i|m|e|f), i_ovvv(m|a|f|e))
                   
                   - contract(e, t1(i|e), 
                            contract(m, t1(m|a), f_ov(m|e)))
                   - contract(m, t1(m|a), 
                            contract(e|n, t1(n|e), i_ooov(m|n|i|e)))
                   + contract(e, t1(i|e), 
                            contract(f|m, t1(m|f), i_ovvv(m|a|f|e)))
                   + contract(f|n, t1(n|f), 
                            contract(m|e, t2(i|m|a|e), i_oovv(m|n|e|f)))
                   - 0.5 * contract(e, t1(i|e), 
                            (contract(f|m|n, t2(m|n|a|f), i_oovv(m|n|e|f))))

                   - 0.5 * contract(m, t1(m|a), 
                            contract(e|f|n, t2(i|n|e|f), i_oovv(m|n|e|f)))
                   - contract(m, t1(m|a), 
                        contract(e, t1(i|e), 
                            contract(f|n, t1(n|f), i_oovv(m|n|e|f))))
                  ;

    return t1_new;

};


/**
 * @brief Action to be applied to CCSD T2 amplitudes.
 * 
 */
container<4,double> action_ccsd_t2 (
    container<2, double> &t1,     ///< CCSD T1 amplitude
    container<4, double> &t2,     ///< CCSD T2 amplitude
    container<2, double> &f_oo,   ///< fock matrix - oo block
    container<2, double> &f_ov,   ///< fock matrix - ov block
    container<2, double> &f_vv,   ///< fock matrix - vv block
    container<4, double> &i_oooo, ///> anti-symmetrized integral <ij||kl>
    container<4, double> &i_ooov, ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv, ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov, ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_ovvv, ///> anti-symmetrized integral <ia||bc>
    container<4, double> &i_vvvv  ///> anti-symmetrized integral <ab||cd>
) {

    container<4, double> t2_new(t2.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;
#if 1
        t2_new(i|j|a|b) = i_oovv(i|j|a|b)
                   + asymm(i, j, contract(e, t1(i|e), i_ovvv(j|e|b|a)))
                   - asymm(a, b, contract(m, t1(m|a), i_ooov(i|j|m|b)))
                   + asymm(a, b, contract(e, t2(i|j|a|e), f_vv(b|e)))
                   - asymm(i, j, contract(m, t2(i|m|a|b), f_oo(j|m)))

                   + 0.5 * contract(e|f, t2(i|j|e|f), i_vvvv(a|b|e|f))
                   - asymm(a, b, asymm(i, j, 
                        contract(e|m, t2(i|m|a|e), i_ovov(m|b|j|e))))
                   + 0.5 * contract(m|n, t2(m|n|a|b), i_oooo(m|n|i|j))
                   + 0.5 * asymm(i, j, 
                        contract(e, t1(i|e), 
                            contract(f, t1(j|f), i_vvvv(a|b|e|f))))

                   - asymm(a, b, asymm(i, j, 
                        contract(e, t1(i|e), 
                            contract(m, t1(m|b), i_ovov(j|e|m|a)))))
                   + 0.5 * asymm(a, b, 
                        contract(m, t1(m|a), 
                            contract(n, t1(n|b), i_oooo(i|j|m|n))))
                   - asymm(a, b, 
                        contract(m, t1(m|a), 
                            contract(e, t2(i|j|e|b), f_ov(m|e))))

                   - asymm(i, j, 
                        contract(e, t1(i|e), 
                            contract(m, t2(m|j|a|b), f_ov(m|e))))
                   - asymm(i, j, 
                        contract(n|e, t1(n|e), 
                            contract(m, t2(m|j|a|b), i_ooov(m|n|i|e))))
                   + asymm(a, b, 
                        contract(e, t2(i|j|e|b),  
                            contract(f|m, t1(m|f), i_ovvv(m|a|f|e))))

                   - 0.5 * asymm(a, b, 
                        contract(m, t1(m|b),  
                            contract(e|f, t2(i|j|e|f), i_ovvv(m|a|f|e))))
                   + 0.5 * asymm(i, j, 
                        contract(e, t1(j|e),  
                            contract(m|n, t2(m|n|a|b), i_ooov(m|n|i|e))))
                   - asymm(a, b, asymm(i, j, 
                        contract(m, t1(m|a), 
                            contract(e|n, t2(n|j|e|b), i_ooov(m|n|i|e)))))

                   + asymm(a, b, asymm(i, j, 
                        contract(e, t1(i|e), 
                            contract(f|m, t2(m|j|f|b), i_ovvv(m|a|f|e)))))
                   - 0.5 * asymm(a, b, asymm(i, j, 
                        contract(e, t1(i|e), 
                            contract(f, t1(j|f), 
                                contract(m, t1(m|b), i_ovvv(m|a|f|e))))))

                   + 0.5 * asymm(a, b, asymm(i, j, 
                        contract(m, t1(m|a), 
                            contract(n, t1(n|b), 
                                contract(e, t1(j|e), i_ooov(m|n|i|e))))))
                   - 0.5 * asymm(i, j, 
                        contract(n, t2(n|j|a|b), 
                            contract(m|e|f, i_oovv(m|n|e|f), t2(i|m|f|e))))

                   - 0.5 * asymm(a, b, 
                        contract(e, t2(i|j|e|b), 
                            contract(f|m|n, t2(m|n|a|f), i_oovv(m|n|e|f))))
                   + 0.5 * asymm(a, b, asymm(i, j, 
                        contract(m|e, t2(i|m|a|e), 
                            contract(n|f, t2(j|n|b|f), i_oovv(m|n|e|f)))))
                   + 0.25 * contract(e|f, t2(i|j|e|f), 
                        contract(m|n, t2(m|n|a|b), i_oovv(m|n|e|f)))


                   - asymm(i, j, 
                        contract(f|n, t1(n|f), 
                            contract(e, t1(i|e), 
                                contract(m, t2(m|j|a|b), i_oovv(m|n|e|f)))))
                   - asymm(a, b, 
                        contract(m, t1(m|a), 
                            contract(f|n, t1(n|f), 
                                contract(e, t2(i|j|e|b), i_oovv(m|n|e|f)))))
                   + 0.25 * asymm(i, j, 
                        contract(e, t1(i|e), 
                            contract(f, t1(j|f), 
                                contract(m|n, t2(m|n|a|b), i_oovv(m|n|e|f)))))

                   + 0.25 * asymm(a, b, 
                        contract(m, t1(m|a), 
                            contract(n, t1(n|b), 
                                contract(e|f, t2(i|j|e|f), i_oovv(m|n|e|f)))))
                   - asymm(a, b, asymm(i, j, 
                        contract(e, t1(i|e), 
                            contract(m, t1(m|a), 
                                contract(f|n, t2(n|j|f|b), i_oovv(m|n|e|f))))))

                   + 0.25 * asymm(a, b, asymm(i, j, 
                        contract(m, t1(m|a), 
                            contract(n, t1(n|b), 
                                contract(e, t1(i|e), 
                                    contract(f, t1(j|f), i_oovv(m|n|e|f)))))))
#endif                                    
                                    ;
    return t2_new;
}        

#endif // ACTION_CCSD_H
