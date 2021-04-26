#ifndef  AMPLITUDES_H
#define  AMPLITUDES_H

#include "supercontainer.h"

enum ampl {t1, t2};

template<typename T=double>
class amplitudes : public supercontainer<T> {
public:
  std::string str(ampl key) {
    std::string str;
    switch (key) {
      case t1: str = "t1"; break;
      case t2: str = "t2"; break;
    }
    return str;
  }

  void set(ampl key, const container<2,T> &c2) { supercontainer<T>::set(str(key), c2); };
  void set(ampl key, const container<4,T> &c4) { supercontainer<T>::set(str(key), c4); };

  container<2,T>& m2get(ampl key) { return supercontainer<T>::m2get(str(key)); };
  container<4,T>& m4get(ampl key) { return supercontainer<T>::m4get(str(key)); };

};

#endif // AMPLITUDES_H
