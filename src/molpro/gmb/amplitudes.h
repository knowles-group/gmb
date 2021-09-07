#ifndef GMB_AMPLITUDES_H
#define GMB_AMPLITUDES_H

#include "supercontainer.h"
#include "utils.h"
#include <numeric>

enum ampl {t1, t2, r1, r2};

template<typename T=double>
class amplitudes : public supercontainer<T> {
private:
  bool triplets{false};
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
  amplitudes& operator=(const std::map<size_t,value_type>& source) {

    size_t count{0}, nsing{0};
    bool doubles{false};
    size_t no{0}, nv{0};
    std::vector<size_t> v_no, v_nv;
    
    for (auto &&im2 : this->m_m2) {
      gmb::zero(*im2.second);
      libtensor::block_tensor_wr_ctrl<2, double> ctrl(*im2.second);
      for (auto &imin : source) {
        libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
        for (libtensor::orbit_list<2, double>::iterator it = ol.begin(); it != ol.end(); it++) {
          libtensor::index<2> bidx;
          ol.get_index(it, bidx);
          if (bidx[0] != bidx[1])
            continue;
          if (bidx[0] == 1)
            count = 0;
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
        if (count < imin.first)
          {
            doubles = true;
            nsing = count;
            
            // get dimensions (#occupied & #virtual)
            libtensor::block_tensor_rd_i<2, double> &bt(*im2.second);
            const libtensor::dimensions<2> &dims = bt.get_bis().get_dims();
            no = dims.get_dim(0);
            nv = dims.get_dim(1);
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
      }
    }

    if (doubles) {

      bool first{true}, ss{false};
      int fact{1};
      size_t ii{0}, jj{0}, aa{0}, bb{0};

    for (auto &&im4 : this->m_m4) {
      gmb::zero(*im4.second);
      libtensor::block_tensor_wr_ctrl<4, double> ctrl(*im4.second);
      for (auto &imin : source) {
        libtensor::orbit_list<4, double> ol(ctrl.req_const_symmetry());

        for (libtensor::orbit_list<4, double>::iterator it = ol.begin(); it != ol.end(); it++) {
          
          libtensor::index<4> bidx;
          ol.get_index(it, bidx);


          if ((bidx[0] != bidx[2] || bidx[1] != bidx[3]) && (bidx[0] != bidx[3] || bidx[1] != bidx[2])
          || ( (bidx[0] > beta) || (bidx[1] > beta) || (bidx[2] > beta) || (bidx[3] > beta)) )
              continue;
          if (ss) {
            if ( (bidx[0] != beta) || (bidx[1] != beta) || (bidx[2] != beta) || (bidx[3] != beta)) {
              continue;
            }
          } else if (!ss) {
            if (bidx[0] == beta && bidx[1] == beta && bidx[2] == beta && bidx[3] == beta) {
              continue;
            } else if (bidx[0] != bidx[2]) {
              fact = -1;
            } else {
              fact = 1;
            } 
          }

          libtensor::dense_tensor_wr_i<4, double> &blk = ctrl.req_block(bidx);
          libtensor::dense_tensor_wr_ctrl<4, double> tc(blk);
          const libtensor::dimensions<4> &tdims = blk.get_dims();
          double *ptr = tc.req_dataptr();
          for (size_t offset = 0; offset < tdims.get_size(); offset++) {
              size_t i = offset / (v_no[bidx[1]]*v_nv[bidx[2]]*v_nv[bidx[3]]);
              size_t j = (offset - i*v_no[bidx[1]]*v_nv[bidx[2]]*v_nv[bidx[3]]) / (v_nv[bidx[2]]*v_nv[bidx[3]]);
              size_t a = (offset - j*v_nv[bidx[2]]*v_nv[bidx[3]] - i*v_no[bidx[1]]*v_nv[bidx[2]]*v_nv[bidx[3]]) / v_nv[bidx[3]];
              size_t b = offset - a*v_nv[bidx[3]] - j*v_nv[bidx[2]]*v_nv[bidx[3]] - i*v_no[bidx[1]]*v_nv[bidx[2]]*v_nv[bidx[3]];
              if (first) {
                if (count == imin.first) {
                  ptr[offset] = imin.second;
                  ii = i;
                  jj = j;
                  aa = a;
                  bb = b;
                  first = false;
                  if (bidx[0] == alpha && bidx[1] == alpha && bidx[2] == alpha && bidx[3] == alpha)
                    ss = true;
                  if (a > b) {
                    ptr[gmb::get_offset(i, j, b, a, v_no[bidx[1]], v_nv[bidx[2]], v_nv[bidx[3]])] = - fact;
                  } else if (i > j){
                    ptr[gmb::get_offset(j, i, a, b, v_no[bidx[1]], v_nv[bidx[2]], v_nv[bidx[3]])] = - fact;
                  }
                }
                ++count;
              } 
              else {
                if (i == ii && j == jj && a == aa && b == bb) {
                  ptr[offset] = fact;
                }
                else if (i == jj && j == ii && a == aa && b == bb){
                  ptr[offset] = - fact;
                }
                else if (i == ii && j == jj && a == bb && b == aa){
                  ptr[offset] = - fact;
                }
                else if (i == jj && j == ii && a == bb && b == aa){
                  ptr[offset] =  fact;
                }
              }
            }
          tc.ret_dataptr(ptr);
          ctrl.ret_block(bidx);
        }
      }
    }
    }
    return *this;
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
    // number of roots needed
    for (size_t ir = 0; ir < n; ir++) 
      m.insert(std::make_pair(1e6+ir, 1e6+ir));

    size_t count{0};

    // singles
    for (auto &&im2 : this->m_m2) {
      libtensor::block_tensor_rd_ctrl<2, double> ctrl(*im2.second);
      libtensor::orbit_list<2, double> ol(ctrl.req_const_symmetry());
      for (libtensor::orbit_list<2, double>::iterator it = ol.begin(); it != ol.end(); it++) {
        libtensor::index<2> bidx;
        ol.get_index(it, bidx);
        if (bidx[0] != bidx[1] || bidx[1] == 1)
          continue;
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
      }
    }

    // doubles
    if (m.rbegin()->first >= 1e6)
      molpro::cout << "Warning: Guess vector will be generated with double excitations; some excited states might not be reasonable.\n";
    for (auto &&im4 : this->m_m4) {
      libtensor::block_tensor_rd_ctrl<4, double> ctrl(*im4.second);
      libtensor::orbit_list<4, double> ol(ctrl.req_const_symmetry());
      for (libtensor::orbit_list<4, double>::iterator it = ol.begin(); it != ol.end(); it++) {
        libtensor::index<4> bidx;
        ol.get_index(it, bidx);
        if ((bidx[0] != bidx[2] || bidx[1] != bidx[3]) && (bidx[0] != bidx[3] || bidx[1] != bidx[2])
          || ( (bidx[0] > beta) || (bidx[1] > beta) || (bidx[2] > beta) || (bidx[3] > beta)) )
              continue;
        libtensor::dense_tensor_rd_i<4, double> &blk = ctrl.req_const_block(bidx);
        libtensor::dense_tensor_rd_ctrl<4, double> tc(blk);
        const libtensor::dimensions<4> &tdims = blk.get_dims();
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
      }
    }
    return m;
  }

};

#endif //GMB_AMPLITUDES_H
