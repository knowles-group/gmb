#ifndef GMB_SRC_MOLPRO_GMB_CCSD_ENERGY_H
#define GMB_SRC_MOLPRO_GMB_CCSD_ENERGY_H

#include <libtensor/libtensor.h>

/**
 * @brief CCSD energy
 * 
 */
double ccsd_energy (
    libtensor::any_tensor<2, double> &t1,     ///< CCSD T1
    libtensor::any_tensor<4, double> &t2,     ///< CCSD T2
    libtensor::any_tensor<2, double> &f_ov,   ///< fock matrix - ov block
    libtensor::any_tensor<4, double> &i_oovv  ///> anti-symmetrized integral <ij||ab>
) {
    
    libtensor::letter i,j,m,n,a,b,e,f;

    double energy = dot_product(t1(m|e), f_ov(m|e))
               + 0.5 * dot_product(t1(m|e)*t1(n|f), i_oovv(m|n|e|f))
               + 0.25 * dot_product(t2(m|n|e|f), i_oovv(m|n|e|f));
    
   return energy;
}

#endif // GMB_SRC_MOLPRO_GMB_CCSD_ENERGY_H
