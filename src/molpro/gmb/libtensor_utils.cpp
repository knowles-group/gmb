#include "libtensor_utils.h"

using namespace libtensor;


namespace gmb {
  template <typename T>
  T dot_prod(libtensor::any_tensor<1,T> &t1,
             libtensor::any_tensor<1,T> &t2) {
    libtensor::letter p;
    T product = dot_product(t2(p), t1(p));
    return product;
  }
  template <typename T>
  T dot_prod(libtensor::any_tensor<2,T> &t1,
             libtensor::any_tensor<2,T> &t2) {
    libtensor::letter p,q;
    T product = dot_product(t2(p|q), t1(p|q));
    return product;
  }
  template <typename T>
  T dot_prod(libtensor::any_tensor<4,T> &t1,
             libtensor::any_tensor<4,T> &t2) {
    libtensor::letter p,q,r,s;
    T product = dot_product(t2(p|q|r|s), t1(p|q|r|s));
    return product;
  }

  template <typename T>
  void copy(libtensor::any_tensor<1,T> &t1,
             libtensor::expr_lhs<1,T> &t2) {
    libtensor::letter p;
    t2(p) = t1(p);
  }

  template <typename T>
  void copy(libtensor::any_tensor<2,T> &t1,
             libtensor::expr_lhs<2,T> &t2) {
    libtensor::letter p, q;
    t2(p|q) = t1(p|q);
  }
                 
  template <typename T>
  void copy(libtensor::any_tensor<4,T> &t1,
            libtensor::expr_lhs<4,T> &t2) {
    libtensor::letter p,q,r,s;
    t2(p|q|r|s) = t1(p|q|r|s);
  }

    template <size_t N, typename T>
        bool comp(const libtensor::btensor<N,T> &tref, 
                  const libtensor::btensor<N,T> &t) {
          libtensor::bto_compare<N,T> cmp(const_cast<libtensor::btensor<N,T>&>(tref), const_cast<libtensor::btensor<N,T>&>(t), 1e-12);
          if(!cmp.compare()) {
              std::ostringstream str;
              cmp.tostr(str);
              std::cout << str.str() << std::endl;
          }
          return cmp.compare();
    }


  template <typename T>
  void compute_axpy(T a,
                    libtensor::any_tensor<1,T> &x,
                    libtensor::expr_lhs<1,T> &y) {
    libtensor::letter p;
    y(p) += a*x(p);
  }
  template <typename T>
  void compute_axpy(T a,
                    libtensor::any_tensor<2,T> &x,
                    libtensor::expr_lhs<2,T> &y) {
    libtensor::letter p,q;
    y(p|q) += a*x(p|q);
  }

  template <typename T>
  void compute_axpy(T a,
                    libtensor::any_tensor<4,T> &x,
                    libtensor::expr_lhs<4,T> &y) {
    libtensor::letter p,q,r,s;
    y(p|q|r|s) += a*x(p|q|r|s);
  }



  void set_sym_pp(libtensor::btensor<2,double> &tensor) {
    libtensor::block_tensor_wr_ctrl<2, double> ctrl(tensor);
    libtensor::symmetry<2, double> &sym = ctrl.req_symmetry();
    // permutational symmetry
    libtensor::permutation<2> p01; p01.permute(0, 1);
    libtensor::scalar_transf<double> tr(1.0);
    libtensor::se_perm<2, double> se(p01, tr);
    sym.insert(se);
  }

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

}

void get_tensor_dimensions(libtensor::btensor<2,double> &tensor, std::vector<size_t> &v_no, std::vector<size_t> &v_nv) {
  // get dimensions (#occupied & #virtual)
  libtensor::block_tensor_rd_i<2, double> &bt(tensor);
  const libtensor::dimensions<2> &dims = bt.get_bis().get_dims();
  auto no = dims.get_dim(0);
  auto nv = dims.get_dim(1);
  auto bis = bt.get_bis();
  // occupied
  const libtensor::split_points &spl_o = bis.get_splits(0);
  for (size_t i = 0; i < spl_o.get_num_points(); i++){
    if (i == 0)
      v_no.push_back(spl_o[i]);
    else 
      v_no.emplace_back(spl_o[i]-spl_o[i-1]);
  }
  v_no.emplace_back(no-std::accumulate(v_no.cbegin(),v_no.cend(),0));
  // virtual
  const libtensor::split_points &spl_v = bis.get_splits(1);
  for (size_t i = 0; i < spl_v.get_num_points(); i++) {
    if (i == 0)
      v_nv.push_back(spl_v[i]);
    else 
      v_nv.emplace_back(spl_v[i]-spl_v[i-1]);
  }
  v_nv.emplace_back(nv-std::accumulate(v_nv.cbegin(),v_nv.cend(),0));
}

template void copy(libtensor::any_tensor<1,double>&,
                   libtensor::expr_lhs<1,double>&); 
template void copy(libtensor::any_tensor<2,double>&,
                   libtensor::expr_lhs<2,double>&); 
template void copy(libtensor::any_tensor<4,double>&,
                   libtensor::expr_lhs<4,double>&); 

template bool comp<2, double> (const libtensor::btensor<2,double> &tref, 
                               const libtensor::btensor<2,double> &t);
template bool comp<4, double> (const libtensor::btensor<4,double> &tref, 
                               const libtensor::btensor<4,double> &t);

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

} // namespace gmb

