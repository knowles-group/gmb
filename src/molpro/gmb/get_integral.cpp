#include "get_integral.h"
#include "utils.h"

bool noel{false}; // no electrons
bool nobeta{false}; // nobeta
bool shift{true}; // shift
bool help{false}; // print stuff
bool zerophoton{false}; // to test - zeroing photonic parts

// // polaritonic parameters  
// extern std::unique_ptr<polariton> ppol;

// get nuclear energy
double get_integral(const std::string &filename) {
  gmb::check_file(filename);
  molpro::FCIdump dump(filename);
  int i, j, k, l;
  double integral(0.0);
  molpro::FCIdump::integralType type;
  dump.rewind();
  while ((type = dump.nextIntegral(i, j, k, l, integral)) != molpro::FCIdump::endOfFile) {
    if (type == molpro::FCIdump::I0)
      if (false) std::cout << "found " <<
                "scalar integral " << integral<< "\n";
  }
  return integral;
}


