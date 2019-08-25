#pragma once

#include <vector>
#include <functional>
#include <numeric>

using VecD = std::vector<double>;

class Scheme;
std::ostream& operator<<(std::ostream& os, const VecD& v);
std::ostream& operator<<(std::ostream& os, const Scheme& s);
VecD operator+(const VecD& lhs, const VecD& rhs);
void operator+=(VecD& lhs, const VecD& rhs);
VecD operator*(const double& lhs, const VecD& rhs);
struct ScalarTimesVector { VecD operator()(const double& lhs, const VecD& rhs); };
struct VectorPlusVector { VecD operator()(VecD const& lhs, VecD const& rhs); };


class Scheme {
  public:
  Scheme() {};
  Scheme(int n, double dt=1e-2);

  int n;
  VecD b;
  VecD c;
  VecD a_lower;
  // std::vector<VecD> a;
  double dt;

  void init();
  void mutate(double prob, double scaling);
  std::pair<VecD, std::vector<VecD>> run(const VecD& x0, std::function<VecD(double,VecD)> f, double t_end) const;
  static Scheme generate(const Scheme&, const Scheme&);

  private:
  std::vector<VecD> calcKVec(double t, const VecD& x, std::function<VecD(double,VecD)> f) const;
  std::vector<VecD> A_full_from_lower(VecD const& a_lower) const;
};

namespace std {
  template<> struct less<Scheme> {
    bool operator() (Scheme const& lhs, Scheme const& rhs) const {
      return std::accumulate(lhs.a_lower.begin(), lhs.a_lower.end(), 0.)
           + std::accumulate(lhs.b.begin(), lhs.b.end(), 0.)
           + std::accumulate(lhs.c.begin(), lhs.c.end(), 0.)
           < std::accumulate(rhs.a_lower.begin(), rhs.a_lower.end(), 0.)
           + std::accumulate(rhs.b.begin(), rhs.b.end(), 0.)
           + std::accumulate(rhs.c.begin(), rhs.c.end(), 0.);
    }
  };
}