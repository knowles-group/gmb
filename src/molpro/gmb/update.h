#ifndef GMB_UPDATE_H
#define GMB_UPDATE_H

#include "container.h"

/**
 * @brief Update T1 amplitudes.
 * 
 */
container<2, double> update_t1 (
    container<2, double> &f_oo,  ///< fock matrix - oo block
    container<2, double> &f_vv,  ///< fock matrix - vv block
    container<2, double> &t1     ///< T1 amplitude
) {

    container<2, double> t1_new(t1.get_space());

    libtensor::letter i,j,k,a,b,c;

    t1_new(i|a) = - div(t1(i|a), dirsum(
                  diag(i, i|j, f_oo(i|j)), -diag(a, a|b, f_vv(a|b))));   
    
    return t1_new;

};

/**
 * @brief Update T2 amplitudes.
 * 
 */
container<4, double> update_t2 (
    container<2, double> &f_oo,  ///< fock matrix - oo block
    container<2, double> &f_vv,  ///< fock matrix - vv block
    container<4, double> &t2      ///< T2 amplitude
) { 

    container<4, double> t2_new(t2.get_space());

    libtensor::letter i,j,k,a,b,c;

    t2_new(i|j|a|b) = - div(t2(i|j|a|b), dirsum(
                    dirsum(diag(i, i|j, f_oo(i|j)), diag(j, j|i, f_oo(j|i))),
                    dirsum(-diag(a, a|b, f_vv(a|b)), -diag(b, b|a, f_vv(b|a)))));   

    return t2_new;
};

#endif //GMB_UPDATE_H
