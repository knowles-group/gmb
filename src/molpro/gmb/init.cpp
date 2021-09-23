#include "init.h"
#include "expressions/fock_xx.h"
#include "expressions/diag_xx.h"

namespace gmb {

  void init(const std::string &filename, const std::string &method, hamiltonian<> &ham, const std::vector<std::unique_ptr<polariton>> &v_ppol) {

    // Two-particle integrals <pq||rs> 
    auto int_oooo = get_i(filename,v_ppol, o, o, o, o);
    auto int_oovv = get_i(filename,v_ppol, o, o, v, v);
    auto int_ovov = get_i(filename,v_ppol, o, v, o, v);
    ham.set(i_oooo, int_oooo);
    ham.set(i_oovv, int_oovv);
    ham.set(i_ovov, int_ovov);

    // One-particle integrals
    auto h1_oo = get_integral(filename,filename,v_ppol,o,o);
    auto h1_vv = get_integral(filename,filename,v_ppol,v,v);

    // get fock matrix
    auto d_oo = diag_xx(h1_oo);
    ham.set(f_oo, fock_xx(d_oo, h1_oo, int_oooo));
    ham.set(f_vv, fock_xx(d_oo, h1_vv, int_ovov));

    // for preconditioner
    auto int_oooo_e = get_i(filename,v_ppol, o, o, o, o, false);
    auto int_ovov_e = get_i(filename,v_ppol, o, v, o, v, false);
    ham.set(f_oo_e, fock_xx(d_oo, h1_oo, int_oooo_e));
    ham.set(f_vv_e, fock_xx(d_oo, h1_vv, int_ovov_e));

    std::ostringstream ss;
    ss << "\nHartree-Fock Orbitals\n\n"
       << "Orbital    Energy (Ha)\n\n";

    readf(ham.m2get(f_oo), ss, 'O');
    readf(ham.m2get(f_vv), ss, 'V');
    molpro::cout << ss.str();


    if (method.find("cc") != std::string::npos) {

      auto int_vvvv = get_i(filename, v_ppol, v, v, v, v);
      ham.set(i_vvvv, int_vvvv);

      if (method.find("ccsd") != std::string::npos) {
        auto h1_ov = get_integral(filename,filename,v_ppol,o,v);
        auto int_ooov = get_i(filename, v_ppol, o, o, o, v);
        auto int_ovvv = get_i(filename, v_ppol, o, v, v, v);

        ham.set(i_vvvv, int_vvvv);
        ham.set(i_ooov, int_ooov);
        ham.set(i_ovvv, int_ovvv);
  
        ham.set(f_ov, fock_xx(d_oo, h1_ov, int_ooov));
      }
    }
  }

  void readf(container<2> &f_xx, std::ostringstream &ss, const char &x) {

      // get dimensions 
      libtensor::block_tensor_rd_i<2, double> &bt(f_xx);
      const libtensor::dimensions<2> &dims = bt.get_bis().get_dims();
      auto norb = dims.get_dim(0);
      auto bis = bt.get_bis();

      std::vector<size_t> v_norb;
      std::vector<double> v_orb;

      const libtensor::split_points &spl = bis.get_splits(0);
      for (size_t i = 0; i < spl.get_num_points(); i++){
        if (i == 0)
          v_norb.push_back(spl[i]);
        else 
          v_norb.emplace_back(spl[i]-spl[i-1]);
      }
      v_norb.emplace_back(norb-std::accumulate(v_norb.cbegin(),v_norb.cend(),0));
      v_orb.reserve(norb);

      // get values
      constexpr size_t N = 2;
      libtensor::block_tensor_rd_ctrl<N, double> ctrl(f_xx);
      libtensor::orbit_list<N, double> ol(ctrl.req_const_symmetry());
      for (libtensor::orbit_list<N, double>::iterator it = ol.begin(); it != ol.end(); it++) {
        libtensor::index<N> bidx;
        ol.get_index(it, bidx);
      if (bidx[0] != bidx[1]) 
        continue;
      libtensor::dense_tensor_rd_i<N, double> &blk = ctrl.req_const_block(bidx);
      libtensor::dense_tensor_rd_ctrl<N, double> tc(blk);
      const libtensor::dimensions<N> &tdims = blk.get_dims();
      const double *ptr = tc.req_const_dataptr();
      for (size_t offset = 0; offset < tdims.get_size(); offset++) {
          size_t i = 1+(offset/v_norb[bidx[1]]);
          size_t j = 1+offset-(offset/v_norb[bidx[1]])*v_norb[bidx[1]];
          if (i == j) {
            ss << std::setw(3) << x << i << gmb::tospin(bidx[0]) << "        "
               << std::setprecision(3) << std::setw(6) <<  ptr[offset] << "\n";
          }
        }
      tc.ret_const_dataptr(ptr);
      ctrl.ret_const_block(bidx);
    }
    ss << "\n";
  }




} // namespace gmb
