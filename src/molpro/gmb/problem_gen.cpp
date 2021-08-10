#include "problem_gen.h"

std::ostream& operator<<(std::ostream& s, const problem_gen& problem) {
  problem.print(s);
  return s;
}