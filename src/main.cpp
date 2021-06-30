#include <molpro/Options.h>
#include <molpro/gmb/gmb.h>
int main(int argc, char *argv[]) {
  molpro::gmb::gmb(molpro::Options("gmb", argc, argv));
  return 0;
}
