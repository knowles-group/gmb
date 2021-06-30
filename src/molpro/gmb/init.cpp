#include "init.h"
#include "expressions/fock_xx.h"
#include "expressions/diag_xx.h"
#include "expressions/anti.h"
#include "expressions/add_d2.h"

// extern std::unique_ptr<polariton> ppol;

namespace gmb {

  void init(const std::string filename, const std::string method, hamiltonian<> &ham, const std::unique_ptr<polariton> &ppol) {

    // One-electron integrals
    auto h1_oo = get_integral(filename,o,o);
    auto h1_vv = get_integral(filename,v,v);

    // Two-electron integrals <pq||rs> 
    auto int_oooo = get_i(filename, o, o, o, o);
    auto int_oovv = get_i(filename, o, o, v, v);
    auto int_ovov = get_i(filename, o, v, o, v);
    
    // add self-energy if needed
    double fact{0};
    std::unique_ptr<container<2>> pd_oo, pd_ov, pd_vo, pd_vv;
    if (ppol != nullptr) { 

      fact = ppol->omega*ppol->gamma*ppol->gamma;

      // add second moment integrals to one-electron part
      auto sm_oo = get_integral(ppol->fname_sm,o,o,false);
      auto sm_vv = get_integral(ppol->fname_sm,v,v,false);
      h1_oo.axpy(fact, sm_oo);
      h1_vv.axpy(fact, sm_vv);

      // add dipole integrals to two-electron part
      pd_oo = std::make_unique<container<2>> (get_integral(ppol->fname_dip,o,o,false));
      pd_ov = std::make_unique<container<2>> (get_integral(ppol->fname_dip,o,v,false));
      pd_vo = std::make_unique<container<2>> (get_integral(ppol->fname_dip,v,o,false));
      pd_vv = std::make_unique<container<2>> (get_integral(ppol->fname_dip,v,v,false));
      add_d2(fact, *pd_oo, *pd_oo, *pd_oo, *pd_oo, int_oooo);
      add_d2(fact, *pd_ov, *pd_ov, *pd_ov, *pd_ov, int_oovv);
      add_d2(fact, *pd_oo, *pd_vv, *pd_vo, *pd_ov, int_ovov);
    }

    ham.set(i_oooo, int_oooo);
    ham.set(i_oovv, int_oovv);
    ham.set(i_ovov, int_ovov);

    // get fock matrix
    auto d_oo = diag_xx(h1_oo);
    ham.set(f_oo, fock_xx(d_oo, h1_oo, int_oooo));
    ham.set(f_vv, fock_xx(d_oo, h1_vv, int_ovov));

    if (method.find("cc") != std::string::npos) {

      auto int_vvvv = get_i(filename, v, v, v, v);
      if (ppol != nullptr) 
        add_d2(fact, *pd_vv, *pd_vv, *pd_vv, *pd_vv, int_vvvv);
      ham.set(i_vvvv, int_vvvv);

      if (method.find("ccsd") != std::string::npos) {
        auto h1_ov = get_integral(filename,o,v);

        auto int_ooov = get_i(filename, o, o, o, v);
        auto int_ovvv = get_i(filename, o, v, v, v);

        if (ppol != nullptr) {
          auto sm_ov = get_integral(ppol->fname_sm,o,v,false);
          h1_ov.axpy(fact, sm_ov);
          add_d2(fact, *pd_oo, *pd_ov, *pd_oo, *pd_ov, int_ooov);
          add_d2(fact, *pd_ov, *pd_vv, *pd_vv, *pd_ov, int_ovvv);
        }
        ham.set(i_ooov, int_ooov);
        ham.set(i_ovvv, int_ovvv);
  
        // getting fock matrix - ov block
        ham.set(f_ov, fock_xx(d_oo, h1_ov, int_ooov));

      }
    }
  }

  container<4,double> get_i(std::string filename, 
                            orb_type o1, 
                            orb_type o2, 
                            orb_type o3, 
                            orb_type o4) {
  
  std::shared_ptr<container<4>> tmp_o1o2o3o4, h2_o1o3o2o4, h2_o1o4o2o3;
  
  h2_o1o3o2o4 = std::make_shared<container<4>> (get_integral(filename, o1, o3, o2, o4)); 

  if (o3 == o4) 
    h2_o1o4o2o3.reset(new container<4>(*h2_o1o3o2o4));
  else 
    h2_o1o4o2o3 = std::make_shared<container<4>> (get_integral(filename, o1, o4, o2, o3)); 
  
  if (o2 == o4 && o2 == o3) 
    tmp_o1o2o3o4.reset(new container<4>(*h2_o1o4o2o3));
  else 
    tmp_o1o2o3o4 = std::make_shared<container<4>> (get_integral(filename, o1, o2, o3, o4)); 

  container<4,double> h2_o1o2o3o4(tmp_o1o2o3o4->get_space());

  // set symmetry
  libtensor::block_tensor_wr_ctrl<4, double> ctrl(h2_o1o2o3o4);
  libtensor::symmetry<4, double> &sym = ctrl.req_symmetry();
  libtensor::scalar_transf<double> tr_sym(1.0);
  libtensor::scalar_transf<double> tr_asym(-1.0);
  if (o1 == o2) {
    libtensor::permutation<4> p01; p01.permute(0, 1);
    libtensor::se_perm<4, double> se_01(p01, tr_asym);
    sym.insert(se_01);
  }
  if (o2 == o3) {
    libtensor::permutation<4> p23; p23.permute(2, 3);
    libtensor::se_perm<4, double> se_23(p23, tr_asym);
    sym.insert(se_23);
  }
  if (o1 == o3 && o2 == o4) {
    libtensor::permutation<4> p0213; p0213.permute(0, 2).permute(1, 3);
    libtensor::se_perm<4, double> se_0213(p0213, tr_sym);
    sym.insert(se_0213);
  }

  anti(h2_o1o2o3o4, *h2_o1o3o2o4, *h2_o1o4o2o3);
  if (false) {
    std::cout << "printing integral " << o1 << o2 << o3 << o4 <<"\n";
    libtensor::bto_print<4, double>(std::cout).perform(h2_o1o2o3o4);
  }  
  return h2_o1o2o3o4;
}



} // namespace gmb

