#ifndef GMB_ID_XXXX_H
#define GMB_ID_XXXX_H

#include "../container.h"

/**
 * @brief Compute identity matrix
 * 
 * Computes identity matrix with same dimensions as a given matrix
 * 
 */

container<4,double> id_xxxx(
    container<4, double> &t_xxxx //> given matrix
) {
    container<4, double> d_xxxx(t_xxxx.get_space()); 
    libtensor::letter p,q,r,s;
    d_xxxx(p|q|r|s) =  set(p|q|r|s, 1.0, set(0.0, d_xxxx(p|q|r|s)));
    return d_xxxx;
}

#endif //GMB_ID_XXXX_H
