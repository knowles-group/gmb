#ifndef GMB_SRC_MOLPRO_GMB_SRC_MOLPRO_GMB_GMB_H_
#define GMB_SRC_MOLPRO_GMB_SRC_MOLPRO_GMB_GMB_H_

#include <molpro/Options.h>
#include <vector>
namespace molpro::gmb {
/*!
 * @brief Carry out gmb calculation
 * @param options
 * @return State energies
 */
std::vector<double>
gmb(const molpro::Options &options = molpro::Options("gmb"));
}
#endif // GMB_SRC_MOLPRO_GMB_GMB_H_
