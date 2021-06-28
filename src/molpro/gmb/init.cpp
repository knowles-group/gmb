#include "init.h"
#include "expressions/fock_xx.h"
#include "expressions/diag_xx.h"
#include "expressions/d2_xxxx.h"
#include "expressions/anti.h"

extern std::unique_ptr<polariton> ppol;
extern std::string filename;

namespace gmb {


  void init(std::string filename, std::string method, hamiltonian<> &ham) {

    // One-electron integrals
    auto h1_oo = get_integral(filename,o,o);
    auto h1_vv = get_integral(filename,v,v);

    // Two-electron integrals <pq||rs> 
    auto int_oooo = get_i(filename, o, o, o, o);
    auto int_oovv = get_i(filename, o, o, v, v);
    auto int_ovov = get_i(filename, o, v, o, v);

    double fact = ppol->omega*ppol->gamma*ppol->gamma;
    auto dip_oo = get_integral(ppol->filename,o,o,false);
    auto dip_ov = get_integral(ppol->filename,o,v,false);
    auto dip_vo = get_integral(ppol->filename,v,o,false);
    auto dip_vv = get_integral(ppol->filename,v,v,false);

    if (ppol != nullptr) {
      libtensor::letter p,q,r,s;
      h1_oo(p|q) += fact*(contract(r,dip_oo(p|r),dip_oo(r|q))
                        + contract(s,dip_ov(p|s),dip_vo(s|q)));
      h1_vv(p|q) += fact*(contract(r,dip_vv(p|r),dip_vv(r|q))
                        + contract(s,dip_vo(p|s),dip_ov(s|q)));
      int_oooo(p|q|r|s) += fact*2.0*(dip_oo(p|r)*dip_oo(q|s)-dip_oo(q|r)*dip_oo(p|s));
      int_oovv(p|q|r|s) += fact*2.0*(dip_ov(p|r)*dip_ov(q|s)-dip_ov(q|r)*dip_ov(p|s));
      int_ovov(p|q|r|s) += fact*2.0*(dip_oo(p|r)*dip_vv(q|s)-dip_vo(q|r)*dip_ov(p|s));
    }

    // getting fock matrix
    auto d_oo = diag_xx(h1_oo);

    ham.set(f_oo, fock_xx(d_oo, h1_oo, int_oooo));
    ham.set(f_vv, fock_xx(d_oo, h1_vv, int_ovov));

    ham.set(i_oooo, int_oooo);
    ham.set(i_oovv, int_oovv);
    ham.set(i_ovov, int_ovov);

    if (method.find("cc") != std::string::npos) {

      auto int_vvvv = get_i(filename, v, v, v, v);

      if (ppol != nullptr) {
        libtensor::letter p,q,r,s;
        int_vvvv(p|q|r|s) += fact*2.0*(dip_vv(p|r)*dip_vv(q|s)-dip_vv(q|r)*dip_vv(p|s));
      }
      ham.set(i_vvvv, int_vvvv);

      if (method.find("ccsd") != std::string::npos) {
        auto h1_ov = get_integral(filename,o,v);

        auto int_ooov = get_i(filename, o, o, o, v);
        auto int_ovvv = get_i(filename, o, v, v, v);

        if (ppol != nullptr) {
          libtensor::letter p,q,r,s;
          h1_ov(p|q) += fact*(contract(r,dip_oo(p|r),dip_ov(r|q))
                            + contract(s,dip_ov(p|s),dip_vv(s|q)));
          int_ooov(p|q|r|s) += fact*2.0*(dip_oo(p|r)*dip_ov(q|s)-dip_oo(q|r)*dip_ov(p|s));
          int_ovvv(p|q|r|s) += fact*2.0*(dip_ov(p|r)*dip_vv(q|s)-dip_vv(q|r)*dip_ov(p|s));
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

