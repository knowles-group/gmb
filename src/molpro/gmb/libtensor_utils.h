#ifndef GMB_LIBTENSOR_UTILS_H
#define GMB_LIBTENSOR_UTILS_H

#include <libtensor/libtensor.h>

// support libtensor functions 
namespace gmb {

  // support libtensor functions for iterative solver

    template <typename T>
    void copy(libtensor::any_tensor<1,T> &t1,
               libtensor::expr_lhs<1,T> &t2);
    template <typename T>
    void copy(libtensor::any_tensor<2,T> &t1,
               libtensor::expr_lhs<2,T> &t2);
                   
    template <typename T>
    void copy(libtensor::any_tensor<4,T> &t1,
               libtensor::expr_lhs<4,T> &t2);

    template <typename T>
    T dot_prod(libtensor::any_tensor<1,T> &t1,
               libtensor::any_tensor<1,T> &t2);

    template <typename T>
    T dot_prod(libtensor::any_tensor<2,T> &t1,
               libtensor::any_tensor<2,T> &t2);

    template <typename T>
    T dot_prod(libtensor::any_tensor<4,T> &t1,
               libtensor::any_tensor<4,T> &t2);

    template <typename T>
    void compute_axpy(T a,
                      libtensor::any_tensor<1,T> &x,
                      libtensor::expr_lhs<1,T> &y);
    template <typename T>
    void compute_axpy(T a,
                      libtensor::any_tensor<2,T> &x,
                      libtensor::expr_lhs<2,T> &y);

    template <typename T>
    void compute_axpy(T a,
                      libtensor::any_tensor<4,T> &x,
                      libtensor::expr_lhs<4,T> &y);

    // set symmetry - probably should be only two functions
    void set_sym_pp(libtensor::btensor<2,double> &tensor);
    void set_sym_pppp(libtensor::btensor<4,double> &tensor);
    void set_sym_pppq(libtensor::btensor<4,double> &tensor);
    void set_sym_ppqq(libtensor::btensor<4,double> &tensor);
    void set_sym_ppqp(libtensor::btensor<4,double> &tensor);
    void set_sym_pqqq(libtensor::btensor<4,double> &tensor);
    void set_sym_pqqp(libtensor::btensor<4,double> &tensor);
    void set_sym_pqpq(libtensor::btensor<4,double> &tensor);
    void set_sym_pqpp(libtensor::btensor<4,double> &tensor);

    // initialize to zero
    void zero(libtensor::btensor<2,double> &tensor);
    void zero(libtensor::btensor<4,double> &tensor);
    
} // namespace gmb
#endif //GMB_LIBTENSOR_UTILS_H
