#ifndef GMB_SRC_MOLPRO_GMB_EOM_CCSD_R0_H
#define GMB_SRC_MOLPRO_GMB_EOM_CCSD_R0_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief Calculate r0
 * 
 */
double eom_ccsd_r0(
    double omega,                  ///< excitation energy
    container<2, double> &r1,      ///< R1 amplitude
    container<4, double> &r2,      ///< R2 amplitude
    container<2, double> &if_ov,   ///< EOM intermediate - ov block
    container<4, double> &i_oovv   ///> anti-symmetrized integral <ij||ab>
) {

    double r0{0};

    libtensor::letter m,n,e,f;

    r0 = ( dot_product(r1(m|e), if_ov(m|e))
       + 0.25*dot_product(r2(m|n|e|f), i_oovv(m|n|e|f)) ) 
       / omega;
                

    return r0;

}

#endif // GMB_SRC_MOLPRO_GMB_EOM_CCSD_R0_H
