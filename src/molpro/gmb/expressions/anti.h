#ifndef GMB_SRC_MOLPRO_GMB_ANTI_H
#define GMB_SRC_MOLPRO_GMB_ANTI_H

#include "../container.h"

/**
 * @brief Computes anti-symmetrized integral
 * 
 */

void anti(
    container<4, double> &i_o1o2o3o4,
    container<4, double> &i_o1o3o2o4,
    container<4, double> &i_o1o4o2o3
) {
    libtensor::letter p,q,r,s;
    // <pq||rs> = <pq|rs> - <pq|sr> = [pr|qs] - [ps|qr]
    i_o1o2o3o4(p|q|r|s) = i_o1o3o2o4(p|r|q|s) - i_o1o4o2o3(p|s|q|r);
}

#endif // GMB_SRC_MOLPRO_GMB_ANTI_H
