#include <gtest/gtest.h>
#include "genetic.h"
#include <cmath>
#include <numeric>
#include <fstream>
#include <fmt/format.h>
#include <fmt/ostream.h>

TEST(SimpsonRule, SinCos) {
  Scheme s(3);
  s.c = {0., 1./2, 1.};
  s.b = {1./6, 4./6, 1./6};
  s.a_lower = {1./2, -1., 2.};
  s.dt = 1e-3;

  auto [timesteps, x] = s.run({0.,1.}, [](double t, const VecD& v) -> VecD { return {v[1], -v[0]}; } , M_PI);

  for (size_t t_i = 0; t_i < timesteps.size(); t_i++) {
    double t = timesteps[t_i];
    ASSERT_NEAR(norm_dist<2>(x[t_i], {sin(t), cos(t)}), 0., std::pow(s.dt, 3)*10);
  }
}

TEST(ExplicitEuler, SinCos) {
  Scheme s(1);
  s.c = {0.};
  s.b = {1.};
  s.a_lower = {};
  s.dt = 0.01;

  auto [timesteps, x] = s.run({0.,1.}, [](double t, const VecD& v) -> VecD { return {v[1], -v[0]}; } , 1);

  for (size_t t_i = 0; t_i < timesteps.size(); t_i++) {
    double t = timesteps[t_i];
    ASSERT_NEAR(norm_dist<2>(x[t_i], {sin(t), cos(t)}), 0., std::pow(s.dt, 1)*10);
  }
}

TEST(HeunScheme, SinCos) {
  Scheme s(2);
  s.c = {0., 1.};
  s.b = {1./2, 1./2};
  s.a_lower = {1.};
  s.dt = 0.01;

  auto [timesteps, x] = s.run({0.,1.}, [](double t, const VecD& v) -> VecD { return {v[1], -v[0]}; } , 1);

  for (size_t t_i = 0; t_i < timesteps.size(); t_i++) {
    double t = timesteps[t_i];
    ASSERT_NEAR(norm_dist<2>(x[t_i], {sin(t), cos(t)}), 0., std::pow(s.dt, 2)*10);
  }
}

