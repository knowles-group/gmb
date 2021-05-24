#include "libtest_utils.h"

namespace gmb_test {

    template <size_t N, typename T>
        bool comp(libtensor::btensor<N,T> tref, 
                  libtensor::btensor<N,T> t) {
          libtensor::bto_compare<N,T> cmp(tref, t, 1e-12);
          if(!cmp.compare()) {
              std::ostringstream str;
              cmp.tostr(str);
              std::cout
                  << "In test : "
                  << str.str() << std::endl;
          }
          return cmp.compare();
    }
 
template bool comp<2, double> (libtensor::btensor<2,double> tref, 
                               libtensor::btensor<2,double> t);
template bool comp<4, double> (libtensor::btensor<4,double> tref, 
                               libtensor::btensor<4,double> t);
       
    
} // namespace gmb_test
