#include <iostream>
#include "genetic.h"

std::ostream& operator<<(std::ostream& os, const VecD& v) {
  os << "[";
  for (double el : v) {
    os << el << " ";
  }
  os << "]";
  return os;
}

std::ostream& operator<<(std::ostream& os, const Scheme& s) {
  os << "Scheme summary:" << std::endl;
  os << "a_lower: " << s.a_lower << std::endl;
  os << "b: " << s.b << std::endl;
  os << "c: " << s.c << std::endl;
  os << "dt: " << s.dt << std::endl;
  return os;
}

VecD operator+(const VecD& lhs, const VecD& rhs) {
  if (lhs.size() != rhs.size()) throw std::runtime_error("Unequal sizes n in operator+!");
  VecD result;
  std::transform(lhs.begin(), lhs.end(), rhs.begin(),
      std::back_inserter(result), std::plus<double>());
  return result;
}

void operator+=(VecD& lhs, const VecD& rhs) {
  if (lhs.size() != rhs.size()) throw std::runtime_error("Unequal sizes n operator+=!");
  lhs = lhs+rhs;
}

VecD operator*(const double& lhs, const VecD& rhs) {
  VecD result;
  std::transform(rhs.begin(), rhs.end(), std::back_inserter(result),
      [lhs](double el)->double {return lhs*el;});
  return result;
}

VecD ScalarTimesVector::operator()(const double& lhs, const VecD& rhs) {
return lhs*rhs;
}

VecD VectorPlusVector::operator()(VecD const& lhs, VecD const& rhs) {
    return lhs+rhs;
}