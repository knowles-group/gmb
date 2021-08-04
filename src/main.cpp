#include <iostream>
#include <molpro/Options.h>
#include <molpro/gmb/gmb.h>
int main(int argc, char *argv[]) {
  auto energies = molpro::gmb::gmb(molpro::Options("gmb", argc, argv));
  std::cout << "Final state energies:";
  for (const auto &e : energies)
    std::cout << " " << e;
  std::cout << std::endl;
  return 0;
}
