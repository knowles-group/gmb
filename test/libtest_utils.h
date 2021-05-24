#include "libtensor/libtensor.h"
#include <libtensor/block_tensor/bto_compare.h>

/**
 * @brief Support functions to test libtensor objects
 * 
 */
namespace gmb_test {

/**
 * @brief Compare 2 tensors
 * 
 * @tparam N tensor rank
 * @tparam T tensor data type
 * @param tref reference tensor
 * @param t tensor to check
 * @return true if they match
 * @return false otherwise
 */
template <size_t N, typename T>
bool comp(libtensor::btensor<N,T> tref, 
          libtensor::btensor<N,T> t);

} // namespace gmb_test

