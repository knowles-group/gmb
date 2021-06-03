#include "init.h"
#include "expressions/fock_xx.h"
#include "expressions/diag_xx.h"

namespace gmb {
  
  void init(std::string filename, std::string method, hamiltonian<> &ham) {

    // getting integrals <pq||rs> 
    auto int_oooo = get_i(filename, o, o, o, o);
    auto int_oovv = get_i(filename, o, o, v, v);
    auto int_ovov = get_i(filename, o, v, o, v);

    // getting fock matrix - oo block
    auto h1_oo = get_integral(filename,o,o);
    auto d_oo = diag_xx(h1_oo);

    // getting fock matrix - vv block
    auto h1_vv = get_integral(filename,v,v);
    ham.set(f_oo, fock_xx(d_oo, h1_oo, int_oooo));
    ham.set(f_vv, fock_xx(d_oo, h1_vv, int_ovov));

    ham.set(i_oooo, int_oooo);
    ham.set(i_oovv, int_oovv);
    ham.set(i_ovov, int_ovov);

    if (method.substr(0,2) == "CC") {

      auto int_vvvv = get_i(filename, v, v, v, v);
      ham.set(i_vvvv, int_vvvv);

      if (method == "CCSD") {
        auto int_ooov = get_i(filename, o, o, o, v);
        auto int_ovvv = get_i(filename, o, v, v, v);
        ham.set(i_ooov, int_ooov);
        ham.set(i_ovvv, int_ovvv);
  
        // getting fock matrix - ov block
        auto h1_ov = get_integral(filename,o,v);
        ham.set(f_ov, fock_xx(d_oo, h1_ov, int_ooov));

      }
    }
  }
  
  void init_ccpol(std::string filename, std::string method, hamiltonian<> &ham) {
    #if 1
    // getting integrals <pq||rs> 
    auto int_oooo = get_i_pol(filename, o, o, o, o);
    auto int_ovov = get_i_pol(filename, o, v, o, v);
    auto int_oovv = get_i_pol(filename, o, o, v, v);
    // int_oooo.scale(10e-10);
    // int_oovv.scale(10e-10);
    // int_ovov.scale(10e-10);
    // getting fock matrix - oo block
    auto h1_oo = get_integral_pol(filename,o,o);
    auto d_oo = diag_xx(h1_oo);

    // getting fock matrix - vv block
    auto h1_vv = get_integral_pol(filename,v,v);

    // h1_oo.scale(0.0);
    // h1_vv.scale(0.0);
    ham.set(f_oo, fock_xx(d_oo, h1_oo, int_oooo));
    ham.set(f_vv, fock_xx(d_oo, h1_vv, int_ovov));

    ham.set(i_oooo, int_oooo);
    ham.set(i_oovv, int_oovv);
    ham.set(i_ovov, int_ovov);

    if (method.substr(0,2) == "CC") {

      auto int_vvvv = get_i_pol(filename, v, v, v, v);
      // int_vvvv.scale(10e-10);
      ham.set(i_vvvv, int_vvvv);

      if (method == "CCSD") {
        auto int_ooov = get_i_pol(filename, o, o, o, v);
        auto int_ovvv = get_i_pol(filename, o, v, v, v);
    
        // int_ooov.scale(10e-10);
        // int_ovvv.scale(10e-10);

        ham.set(i_ooov, int_ooov);
        ham.set(i_ovvv, int_ovvv);
  
        // getting fock matrix - ov block
        auto h1_ov = get_integral_pol(filename,o,v);
        // h1_ov.scale(0.0);
        ham.set(f_ov, fock_xx(d_oo, h1_ov, int_ooov));

      }
    }
    #endif
  }

} // namespace gmb

