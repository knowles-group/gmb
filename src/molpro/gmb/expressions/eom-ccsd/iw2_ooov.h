#ifndef GMB_SRC_MOLPRO_GMB_EOM_CCSD_IW2_OOOV_H
#define GMB_SRC_MOLPRO_GMB_EOM_CCSD_IW2_OOOV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ooov intermediate
 * 
 */
container<4,double> eom_ccsd_iw2_ooov(
    container<2, double> &t1,       ///< CCSD T1
    container<4, double> &i_ooov,   ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv    ///> anti-symmetrized integral <ij||ab>

) {

    container<4, double> iw_ooov(i_ooov.get_space());

    libtensor::letter a,i,j,k,e;

    iw_ooov(j|k|i|a) = i_ooov(j|k|i|a)
                     + contract(e, t1(i|e), i_oovv(j|k|e|a));
    return iw_ooov;

}

#endif // GMB_SRC_MOLPRO_GMB_EOM_CCSD_IW2_OOOV_H
