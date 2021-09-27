#ifndef GMB_SRC_MOLPRO_GMB_EOM_CCSD_TAU_H
#define GMB_SRC_MOLPRO_GMB_EOM_CCSD_TAU_H

#include <libtensor/libtensor.h>
#include "../../container.h"


/**
 * @brief tau intermediate
 * 
 */
container<4,double> eom_ccsd_tau(
    container<2, double> &t1,     ///< CCSD T1
    container<4, double> &t2      ///< CCSD T2
) {

    container<4, double> tau(t2.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    tau(i|j|a|b) = t2(i|j|a|b) + asymm(i, j, t1(i|a)*t1(j|b));
        // tau(i|j|a|b) =
        //    0.25 * asymm(i,j, asymm(a,b, 
        //       t2(i|j|a|b)
        //     + asymm(a, b, asymm(i, j, 0.5 * t1(i|a) * t1(j|b)))));

    return tau;

}


#endif // GMB_SRC_MOLPRO_GMB_EOM_CCSD_TAU_H
