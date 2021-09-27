#ifndef GMB_SRC_MOLPRO_GMB_ID_XX_H
#define GMB_SRC_MOLPRO_GMB_ID_XX_H

#include "../container.h"

/**
 * @brief Compute identity matrix
 * 
 * Computes identity matrix with same dimensions as a given matrix
 * 
 */

container<2,double> id_xx(
    container<2, double> &t_xx //> given matrix
) {
    container<2, double> d_xx(t_xx.get_space()); 
    libtensor::letter p,q;
    d_xx(p|q) =  set(p|q, 1.0, set(0.0, d_xx(p|q)));
    return d_xx;
}

#endif // GMB_SRC_MOLPRO_GMB_ID_XX_H
