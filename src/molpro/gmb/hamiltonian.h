#ifndef GMB_HAMILTONIAN_H
#define GMB_HAMILTONIAN_H

#include "supercontainer.h"

// one-electron (oe) part
enum oe_part {f_oo, f_ov, f_vv, f_oo_e, f_vv_e};
// two-electron (te) part
enum te_part {i_oooo, i_ooov, i_oovv, i_ovvv, i_ovov, i_vvvv};

template<typename T=double>
class hamiltonian : public supercontainer<T> {
public:
  std::string str(oe_part key) {
    std::string str;
    switch (key) {
      case f_oo: str = "f_oo"; break;
      case f_oo_e: str = "f_oo_e"; break;
      case f_ov: str = "f_ov"; break;
      case f_vv: str = "f_vv"; break;
      case f_vv_e: str = "f_vv_e"; break;
    }
    return str;
  }
  
  std::string str(te_part key) {
    std::string str;
    switch (key) {
      case i_oooo: str = "i_oooo"; break;
      case i_ooov: str = "i_ooov"; break;
      case i_oovv: str = "i_oovv"; break;
      case i_ovov: str = "i_ovov"; break;
      case i_ovvv: str = "i_ovvv"; break;
      case i_vvvv: str = "i_vvvv"; break;
    }
    return str;
  }

  void set(oe_part key, const container<2,T> &c2) { supercontainer<T>::set(str(key), c2); };
  void set(te_part key, const container<4,T> &c4) { supercontainer<T>::set(str(key), c4); };

  container<2,T>& m2get(oe_part key) { return supercontainer<T>::m2get(str(key)); };
  container<4,T>& m4get(te_part key) { return supercontainer<T>::m4get(str(key)); };

};

#endif //GMB_HAMILTONIAN_H
