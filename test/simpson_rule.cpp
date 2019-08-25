#include <gtest/gtest.h>
#include "genetic.h"
#include <cmath>
#include <numeric>
#include <fstream>
#include <fmt/format.h>
#include <fmt/ostream.h>

template<size_t p>
double norm_dist(VecD const& lhs, VecD const& rhs) {
  VecD squared_errors(std::min(lhs.size(), rhs.size()));
  std::transform(lhs.begin(), lhs.end(),
                 rhs.begin(),
                 squared_errors.begin(),
                 [](double lhs, double rhs) {return std::pow(lhs-rhs, p);});
  return std::pow(std::accumulate(squared_errors.begin(), squared_errors.end(), 0.), 1./p);
}

TEST(SimpsonRule, SinCos) {
  Scheme s(3);
  s.c = {0., 1./2, 1.};
  s.b = {1./6, 4./6, 1./6};
  s.a_lower = {1./2, -1., 2.};
  s.dt = 1e-2;

  auto x = s.run({0.,1.}, [](double t, const VecD& v) -> VecD { return {v[1], -v[0]}; } , M_PI);
  std::ofstream sincos("sincos.csv");
  ASSERT_TRUE(sincos.is_open());
  for (auto x_ : x) {
    sincos << x_[0] << "," << x_[1] << std::endl;
    //fmt::print(sincos, "{1},{2}", x_[0], x_[1]);
  }
  sincos.close();

  size_t N_t = (M_PI-0)/s.dt + 1;
  VecD timesteps(N_t); std::generate(timesteps.begin(), timesteps.end(), [n=0, N_t]() mutable {return n++/(N_t-1);});

  for (size_t t_i = 0; t_i < N_t; t_i++) {
    double t = timesteps[t_i];
    ASSERT_NEAR(norm_dist<2>(x[t_i], {sin(t), cos(t)}), 0., s.dt);
  }
  // ASSERT_NEAR(std::sqrt(std::pow(x.back()[0]-std::cos(1),2)+std::pow(x.back()[1]+std::sin(1),2)), 0, 2e-9);
}

TEST(ExplicitEuler, SinCos) {
  Scheme s(1);
  s.c = {0.};
  s.b = {1.};
  s.a_lower = {};
  s.dt = 0.01;

  auto x = s.run({0.,1.}, [](double t, const VecD& v) -> VecD { return {v[1], -v[0]}; } , 1);

  ASSERT_NEAR(std::sqrt(std::pow(x.back()[0]-std::cos(1),2)+std::pow(x.back()[1]+std::sin(1),2)), 0, 2e-2);
}

TEST(HeunScheme, SinCos) {
  Scheme s(2);
  s.c = {0., 1.};
  s.b = {1./2, 1./2};
  s.a_lower = {1.};
  s.dt = 0.01;

  auto x = s.run({0.,1.}, [](double t, const VecD& v) -> VecD { return {v[1], -v[0]}; } , 1);

  ASSERT_NEAR(std::sqrt(std::pow(x.back()[0]-std::cos(1),2)+std::pow(x.back()[1]+std::sin(1),2)), 0, 2e-4);
}

