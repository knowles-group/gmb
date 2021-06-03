#ifndef GMB_AMPLITUDES_H
#define GMB_AMPLITUDES_H

#include "supercontainer.h"

extern std::string filename;
enum ampl {t1, t2, r1, r2};


template<typename T=double>
class amplitudes : public supercontainer<T> {
private:
  bool triplets = false;
public:
  std::string str(ampl key) {
    std::string str;
    switch (key) {
      case t1: str = "t1"; break;
      case t2: str = "t2"; break;
      case r1: str = "r1"; break;
      case r2: str = "r2"; break;
    }
    return str;
  }
  
  using value_type = T;

  using supercontainer<T>::supercontainer;

  amplitudes() 
  : supercontainer<T>() {}

/**
 * @brief Construct a new amplitudes object to be used as guess vector.
 * 
 * The object is constructed from a map which contains what position should be filled with which value.
 * 
 * @param source map containing position and value to be added
 */
  amplitudes(const std::map<size_t,value_type>& source) {
    for (auto &imin : source) {
      auto r1_guess = get_integral(filename,o,v);
      gmb::zero(r1_guess);
      libtensor::block_tensor_wr_ctrl<2, double> ctrl(r1_guess);
      libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
      size_t count(0);
      for (libtensor::orbit_list<2, double>::iterator it = ol.begin();
           it != ol.end(); it++) {
        libtensor::index<2> bidx;
        ol.get_index(it, bidx);
        libtensor::dense_tensor_wr_i<2, double> &blk = ctrl.req_block(bidx);
        libtensor::dense_tensor_wr_ctrl<2, double> tc(blk);
        const libtensor::dimensions<2> &tdims = blk.get_dims();
        double *ptr = tc.req_dataptr();
        for (size_t i = 0; i < tdims.get_size(); i++) {
          if (count == imin.first) 
            ptr[i] = imin.second;
          ++count;
        }
        tc.ret_dataptr(ptr);
        ctrl.ret_block(bidx);
      }
      this->m_m2.insert(std::make_pair("r1", new container<2,T> (r1_guess)));
      auto r2_guess = get_integral(filename,o,o,v,v);
      gmb::zero(r2_guess);
      this->m_m4.insert(std::make_pair("r2", new container<4,T> (r2_guess)));
    }
  }

  void set(ampl key, const container<2,T> &c2) { supercontainer<T>::set(str(key), c2); };
  void set(ampl key, const container<4,T> &c4) { supercontainer<T>::set(str(key), c4); };

  container<2,T>& m2get(ampl key) { return supercontainer<T>::m2get(str(key)); };
  container<4,T>& m4get(ampl key) { return supercontainer<T>::m4get(str(key)); };
  
  /**
   * @brief Selects smallest n values of container. 
   * Returns a map with the value and the position.
   * 
   * @param n number of smallest values needed
   * @param max 
   * @param ignore_sign 
   * @return std::map<size_t, T> 
   */
  std::map<size_t, T> select(size_t n, bool max = false, bool ignore_sign = false) const {
    std::map<size_t, T> m;
    size_t count{0};

    // for now only singles
    for (auto &&im2 : this->m_m2) {
      // number of roots needed
      for (size_t ir = 0; ir < n; ir++) 
        m.insert(std::make_pair(1e6+ir, 1e6+ir));
      libtensor::block_tensor_rd_ctrl<2, double> ctrl(*im2.second);
      libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
      for (libtensor::orbit_list<2, double>::iterator it = ol.begin();
           it != ol.end(); it++) {
        libtensor::index<2> bidx;
        ol.get_index(it, bidx);
        libtensor::dense_tensor_rd_i<2, double> &blk = ctrl.req_const_block(bidx);
        libtensor::dense_tensor_rd_ctrl<2, double> tc(blk);
        const libtensor::dimensions<2> &tdims = blk.get_dims();
        const double *ptr = tc.req_const_dataptr();
        for (size_t i = 0; i < tdims.get_size(); i++) {
            typename std::map<size_t,T>::iterator mmax
              = std::max_element(m.begin(),m.end(),[](const std::pair<size_t,T>& a, const std::pair<size_t,T>& b)->bool{ return a.second < b.second; } );
          if (ptr[i] < mmax->second) {
            m.erase(mmax);
            m.insert(std::make_pair(count, ptr[i]));
          }
          ++count;
        }
        tc.ret_const_dataptr(ptr);
        ctrl.ret_const_block(bidx);
        if (!triplets) break; // only alpha-alpha excitations
      }
    }

    return m;
  }
};

#endif //GMB_AMPLITUDES_H
