#include "libtensor_utils.h"

using namespace libtensor;


namespace bbo {
  template <typename T>
  T dot_prod(libtensor::any_tensor<1,T> &t1,
             libtensor::any_tensor<1,T> &t2) {
    libtensor::letter p;
    T product = dot_product(t2(p), t1(p));
    return product;
  };
  template <typename T>
  T dot_prod(libtensor::any_tensor<2,T> &t1,
             libtensor::any_tensor<2,T> &t2) {
    libtensor::letter p,q;
    T product = dot_product(t2(p|q), t1(p|q));
    return product;
  };
  template <typename T>
  T dot_prod(libtensor::any_tensor<4,T> &t1,
             libtensor::any_tensor<4,T> &t2) {
    libtensor::letter p,q,r,s;
    T product = dot_product(t2(p|q|r|s), t1(p|q|r|s));
    return product;
  };

  template <typename T>
  void copy(libtensor::any_tensor<1,T> &t1,
             libtensor::expr_lhs<1,T> &t2) {
    libtensor::letter p;
    t2(p) = t1(p);
  };
  template <typename T>
  void copy(libtensor::any_tensor<2,T> &t1,
             libtensor::expr_lhs<2,T> &t2) {
    libtensor::letter p, q;
    t2(p|q) = t1(p|q);
  };
                 
  template <typename T>
  void copy(libtensor::any_tensor<4,T> &t1,
            libtensor::expr_lhs<4,T> &t2) {
    libtensor::letter p,q,r,s;
    t2(p|q|r|s) = t1(p|q|r|s);
  };
  template <typename T>
  void compute_axpy(T a,
                    libtensor::any_tensor<1,T> &x,
                    libtensor::expr_lhs<1,T> &y) {
    libtensor::letter p;
    y(p) += a*x(p);
  };
  template <typename T>
  void compute_axpy(T a,
                    libtensor::any_tensor<2,T> &x,
                    libtensor::expr_lhs<2,T> &y) {
    libtensor::letter p,q;
    y(p|q) += a*x(p|q);
  };

  template <typename T>
  void compute_axpy(T a,
                    libtensor::any_tensor<4,T> &x,
                    libtensor::expr_lhs<4,T> &y) {
    libtensor::letter p,q,r,s;
    y(p|q|r|s) += a*x(p|q|r|s);
  };

  void set_sym_pp(libtensor::btensor<2,double> &tensor) {
    libtensor::block_tensor_wr_ctrl<2, double> ctrl(tensor);
    libtensor::symmetry<2, double> &sym = ctrl.req_symmetry();
    // permutational symmetry
    libtensor::permutation<2> p01; p01.permute(0, 1);
    libtensor::scalar_transf<double> tr(1.0);
    libtensor::se_perm<2, double> se(p01, tr);
    sym.insert(se);
  };

  void set_sym_pppp(libtensor::btensor<4,double> &tensor) {
    // Request a control object
    libtensor::block_tensor_wr_ctrl<4, double> ctrl(tensor);
    libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();

    // permutational symmetry
    libtensor::scalar_transf<double> tr(1.0);
    libtensor::permutation<4> p01; p01.permute(0, 1);
    libtensor::se_perm<4, double> se_01(p01, tr);
    sym.insert(se_01);
    libtensor::permutation<4> p23; p23.permute(2, 3);
    libtensor::se_perm<4, double> se_23(p23, tr);
    sym.insert(se_23);
    libtensor::permutation<4> p0213; p0213.permute(0, 2).permute(1, 3);
    libtensor::se_perm<4, double> se_0213(p0213, tr);
    sym.insert(se_0213);
  }

  void set_sym_pppq(libtensor::btensor<4,double> &tensor) {
    // Request a control object
    libtensor::block_tensor_wr_ctrl<4, double> ctrl(tensor);
    libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();

    // permutational symmetry
    libtensor::permutation<4> p01; p01.permute(0, 1);
    libtensor::scalar_transf<double> tr(1.0);
    libtensor::se_perm<4, double> se_01(p01, tr);
    sym.insert(se_01);
  }

  void set_sym_ppqq(libtensor::btensor<4,double> &tensor) {
    // Request a control object
    libtensor::block_tensor_wr_ctrl<4, double> ctrl(tensor);
    libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();

    // permutational symmetry
    libtensor::permutation<4> p01; p01.permute(0, 1);
    libtensor::permutation<4> p23; p23.permute(2, 3);
    libtensor::scalar_transf<double> tr(1.0);
    libtensor::se_perm<4, double> se_01(p01, tr);
    libtensor::se_perm<4, double> se_23(p23, tr);
    sym.insert(se_01);
    sym.insert(se_23);
  }

  void set_sym_ppqp(libtensor::btensor<4,double> &tensor) {
    // Request a control object
    libtensor::block_tensor_wr_ctrl<4, double> ctrl(tensor);
    libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();

    // permutational symmetry
    libtensor::permutation<4> p01; p01.permute(0, 1);
    libtensor::scalar_transf<double> tr(1.0);
    libtensor::se_perm<4, double> se_01(p01, tr);
    sym.insert(se_01);
  }
  void set_sym_pqqq(libtensor::btensor<4,double> &tensor) {
    // Request a control object
    libtensor::block_tensor_wr_ctrl<4, double> ctrl(tensor);
    libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();

    // permutational symmetry
    libtensor::permutation<4> p23; p23.permute(2, 3);
    libtensor::scalar_transf<double> tr(1.0);
    libtensor::se_perm<4, double> se_23(p23, tr);
    sym.insert(se_23);
  }


  void set_sym_pqpq(libtensor::btensor<4,double> &tensor) {
    // Request a control object
    libtensor::block_tensor_wr_ctrl<4, double> ctrl(tensor);
    libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();

    // permutational symmetry

    libtensor::scalar_transf<double> tr(1.0);
    libtensor::permutation<4> p0213; p0213.permute(0, 2).permute(1, 3);
    libtensor::se_perm<4, double> se_0213(p0213, tr);
    sym.insert(se_0213);
  }

  void set_sym_pqpp(libtensor::btensor<4,double> &tensor) {
    // Request a control object
    libtensor::block_tensor_wr_ctrl<4, double> ctrl(tensor);
    libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();

    // permutational symmetry
    libtensor::permutation<4> p23; p23.permute(2, 3);
    libtensor::scalar_transf<double> tr(1.0);
    libtensor::se_perm<4, double> se_23(p23, tr);
    sym.insert(se_23);
  }

  void zero(libtensor::btensor<2,double> &tensor) {
	if (false) std::cout << "copying a two-dimensional tensor\n";

  libtensor::block_tensor_wr_ctrl<2, double> ctrl(tensor);

  // Loop over all canonical blocks using orbit_list
  libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
  // libtensor::orbit_list<4, double> ol_h2(ctrl_h2.req_const_symmetry());
  for (libtensor::orbit_list<2, double>::iterator it = ol.begin();
           it != ol.end(); it++) {
    
    // Obtain the index of the current canonical block
    libtensor::index<2> bidx;
    ol.get_index(it, bidx);
    // Request tensor block from control object
    libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
    // Fill with data 
    libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);

    // Obtain dimensions of tensor block
    const libtensor::dimensions<2> &tdims = blk.get_dims();
    // Request data pointer
    double *ptr = tc.req_dataptr();
    for (size_t i = 0; i < tdims.get_size(); i++) {
        ptr[i] = 0.0;
    }
    // Return data pointer
    tc.ret_dataptr(ptr);

    // Return the tensor block (mark as done)
    ctrl.ret_block(bidx);
  }

  if (false) {
  std::cout << "printing t\n";
  libtensor::bto_print<2, double>(std::cout).perform(tensor);
  }
}

  // initialize to zero
  void zero(libtensor::btensor<4,double> &tensor) {

  libtensor::block_tensor_wr_ctrl<4, double> ctrl(tensor);
  // Loop over all canonical blocks using orbit_list
  libtensor::orbit_list<4, double> ol(ctrl.req_const_symmetry());
  for (libtensor::orbit_list<4, double>::iterator it = ol.begin();
           it != ol.end(); it++) {
    
    // Obtain the index of the current canonical block
    libtensor::index<4> bidx;
    ol.get_index(it, bidx);
    // Request tensor block from control object
    libtensor::dense_tensor_wr_i<4, double> &blk = ctrl.req_block(bidx);
    // Fill with data 
    libtensor::dense_tensor_wr_ctrl<4, double> tc(blk);
    // Obtain dimensions of tensor block
    const libtensor::dimensions<4> &tdims = blk.get_dims();
    // Request data pointer
    double *ptr = tc.req_dataptr();
    for (size_t i = 0; i < tdims.get_size(); i++) {
            ptr[i] = 0.0;    
    }
     
    // Return data pointer
    tc.ret_dataptr(ptr);
    // Return the tensor block (mark as done)
    ctrl.ret_block(bidx);
  }
  if (false) {
  std::cout << "printing t\n";
  libtensor::bto_print<4, double>(std::cout).perform(tensor);
  }
}

template void copy(libtensor::any_tensor<1,double>&,
                   libtensor::expr_lhs<1,double>&); 
template void copy(libtensor::any_tensor<2,double>&,
                   libtensor::expr_lhs<2,double>&); 
template void copy(libtensor::any_tensor<4,double>&,
                   libtensor::expr_lhs<4,double>&); 
template double dot_prod(libtensor::any_tensor<1,double>&,
                         libtensor::any_tensor<1,double>&); 
template double dot_prod(libtensor::any_tensor<2,double>&,
                         libtensor::any_tensor<2,double>&); 
template double dot_prod(libtensor::any_tensor<4,double>&,
                         libtensor::any_tensor<4,double>&);
                         
template void compute_axpy(double,
                           libtensor::any_tensor<1,double>&,
                           libtensor::expr_lhs<1,double>&);
template void compute_axpy(double,
                           libtensor::any_tensor<2,double>&,
                           libtensor::expr_lhs<2,double>&);
template void compute_axpy(double,
                           libtensor::any_tensor<4,double>&,
                           libtensor::expr_lhs<4,double>&);

} // namespace bbo

