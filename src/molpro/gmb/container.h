#ifndef GMB_CONTAINER_H
#define GMB_CONTAINER_H


#include <molpro/iostream.h>
#include <libtensor/libtensor.h>
#include <libtensor/block_tensor/bto_set.h>
#include "libtensor_utils.h"

// using namespace libtensor; 
/**
 * @brief Container class for btensor objects
 * 
 * @tparam N tensor rank
 * @tparam T tensor data type
 */

template <size_t N, typename T=double>
class container : public libtensor::btensor<N,T>,
     virtual public libtensor::block_tensor<N, T, libtensor::allocator> {
private:
  libtensor::bispace<N> m_sp; ///< tensor space
public:

  explicit container(const libtensor::bispace<N> &sp) 
  : m_sp(sp), 
    libtensor::btensor<N,T>(sp), 
    libtensor::block_tensor< N, T, libtensor::allocator >(sp.get_bis())
    {}

  container(const container &c) 
  : m_sp(c.get_space()), 
    libtensor::btensor<N,T>(c.get_space()), 
    libtensor::block_tensor< N, T, libtensor::allocator > (c.get_space().get_bis()) {
    gmb::copy(const_cast<container&> (c),*this);
  }

  container& operator=(const container &rhs) {
    gmb::copy(const_cast<container&>(rhs), *this);
    return *this;
  }

    using value_type = T;

  friend bool operator==(const container<N,T>& lhs, const container<N,T>& rhs) {
    // return gmb::comp(const_cast<container&>(lhs), const_cast<container&>(rhs));
    return gmb::comp(lhs, rhs);
  }


  /**
   * @brief Get the space object
   * 
   * @return const libtensor::bispace<N>& 
   */
  const libtensor::bispace<N>& get_space() const { return m_sp; };

  /**
   * @brief Fill a tensor with a constant value a
   * 
   * @param a constant value to fill the tensor with
   */
  void fill(T a) { libtensor::bto_set<N,T>(a).perform(*this); };

  /**
   * @brief Scale the tensor with constant value a
   * 
   * @param a constant to be scale the tensor by
   */
  void scal(T a) { (*this).scale(a); };

  /**
   * @brief Dot tensor with another tensor
   * 
   * @param x tensor to be dotted with
   * @return T dot product between tensors
   */
  T dot(const container &x) const { return gmb::dot_prod(const_cast<container&> (*this), const_cast<container&> (x)); };

  /**
   * @brief Assigns new value to the tensor Y = a*X + Y
   * 
   * @param a constant value
   * @param x another tensor
   */
  void axpy(T a, const container &x) { gmb::compute_axpy(a, const_cast<container&> (x), (*this)); };
  // void axpy(T a, container x) { gmb::compute_axpy(a, x, (*this)); };

  /**
   * @brief print container
   * 
   */
  void print() { libtensor::bto_print<N, T>(molpro::cout).perform(*this);};


  std::map<size_t, value_type> select_max_dot(size_t n, const container& y) const {throw std::logic_error("container::select_max_dot unimplemented");}



};

#endif //GMB_CONTAINER_H
