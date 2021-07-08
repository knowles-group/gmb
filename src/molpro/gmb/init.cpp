#include "init.h"
#include "expressions/fock_xx.h"
#include "expressions/diag_xx.h"

namespace gmb {

  void init(const std::string &filename, const std::string &method, hamiltonian<> &ham, const std::vector<std::shared_ptr<polariton>> &v_ppol) {

    // Two-particle integrals <pq||rs> 
    auto int_oooo = get_i(filename, v_ppol, o, o, o, o);
    auto int_oovv = get_i(filename, v_ppol, o, o, v, v);
    auto int_ovov = get_i(filename, v_ppol, o, v, o, v);

    ham.set(i_oooo, int_oooo);
    ham.set(i_oovv, int_oovv);
    ham.set(i_ovov, int_ovov);

    // One-particle integrals
    auto h1_oo = get_integral(filename, filename,v_ppol,o,o);
    auto h1_vv = get_integral(filename, filename,v_ppol,v,v);

    // get fock matrix
    auto d_oo = diag_xx(h1_oo);
    ham.set(f_oo, fock_xx(d_oo, h1_oo, int_oooo));
    ham.set(f_vv, fock_xx(d_oo, h1_vv, int_ovov));

    if (method.find("cc") != std::string::npos) {

      auto int_vvvv = get_i(filename, v_ppol, v, v, v, v);
      ham.set(i_vvvv, int_vvvv);

      if (method.find("ccsd") != std::string::npos) {
        auto h1_ov = get_integral(filename, filename,v_ppol,o,v);

        auto int_ooov = get_i(filename, v_ppol, o, o, o, v);
        auto int_ovvv = get_i(filename, v_ppol, o, v, v, v);

        ham.set(i_ooov, int_ooov);
        ham.set(i_ovvv, int_ovvv);
  
        ham.set(f_ov, fock_xx(d_oo, h1_ov, int_ooov));

      }
    }
  }

} // namespace gmb

