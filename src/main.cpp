#include "genetic.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <random>
#include <numeric>
#include <iomanip>

using VecS = std::vector<Scheme>;

std::tuple<double,VecD> evalSinCos(const Scheme& s) {
  double t_end = 1;
  auto result = s.run({1.,0.},[](double t, const VecD& v)->VecD {return {v[1],-v[0]};}, t_end);
  double error = 0;
  int n_t = std::floor((t_end-0)/s.dt);
  double t = 0;
  for (int i = 0; i < n_t; i++) {
    error += std::sqrt(std::pow(result[i][0]-std::cos(t),2)+std::pow(result[i][1]-(-std::sin(t)),2));
    t += s.dt;
  }

  return {error, result[result.size()-1]};
}

std::tuple<double,VecD> eval2DSecondOrder(const Scheme& s) {
  double t_end = 0.1;
  auto result = s.run({1.,1.},[](double t, const VecD& v)->VecD {return {v[1],-10*v[1]+std::pow(t,2)};}, t_end);
  double error = 0;

  int n_t = std::floor((t_end-0)/s.dt)+1;
  double t = 0;
  for (int i = 0; i < n_t; i++) {
    double C1 = -499./5000;
    double C2 = 1-C1;
    double y = std::pow(t,3)/30. - std::pow(t,2)/.100 + t/500. + C1*std::exp(-10.*t) + C2;
    error += std::pow(result[i][0]-y,2);
    t += s.dt;
  }

  return {error, result[result.size()-1]};
}

std::tuple<double,VecD> eval(const Scheme& s) {
  double t_end = 1;
  auto result = s.run({1.,0.},[](double t, const VecD& v)->VecD {return {v[1]*t, -v[0]*t};}, t_end);
  double error = 0;

  int n_t = std::floor((t_end-0)/s.dt)+1;
  double t = 0;
  for (int i = 0; i < n_t; i++) {
    VecD y = {std::cos(t*t/2), -std::sin(t*t/2)};
    error += std::sqrt(std::pow(result[i][0]-y[0],2)+std::pow(result[i][1]-y[1],2));
    t += s.dt;
  }

  return {error, result[result.size()-1]};
}


int main() {
  std::random_device rd{};
  std::mt19937 gen(rd());
  const int sample_size = 100;
  const int epochs = 50;
  const int fittest = 15;
  VecS samples(sample_size,Scheme(3,1e-3));
  std::for_each(samples.begin(), samples.end(), 
      [](Scheme& s) {s.init();});

  // Set one scheme to converged scheme for debug
  /* samples[50].a_lower = {1./2,-1.,2.}; */
  /* samples[50].b = {1./6, 4./6, 1./6}; */
  /* samples[50].c = {0., 1./2, 1.}; */

  // do epochs
  for (int i = 0; i < epochs; i++) {
    // sort by fitness and log
    std::clog << "Sorting samples" << std::endl;
    std::sort(samples.begin(), samples.end(), 
        [](const Scheme& s1, const Scheme& s2){
          auto [fitness1, x1] = evalSinCos(s1);
          auto [fitness2, x2] = evalSinCos(s2); 
          return fitness1 < fitness2;
        });
    std::clog << "Done." << std::endl;
    std::clog << "Evaluating best sample..." << std::endl;
    auto [fitness, x_end] = evalSinCos(samples[0]);
    std::clog << std::setprecision(15);
    std::clog << "Best Score was: " << fitness << "!" << std::endl;
    std::clog << "x = " << x_end << std::endl;
    std::clog << std::setprecision(5);
    std::clog << samples[0] << std::endl;
    std::clog << "#####################" << std::endl << std::endl;;

    // replace unfitt by offspring of fittest population
    for (auto it = std::next(samples.begin(),fittest); it != samples.end(); it++) {
      std::vector<Scheme> parents;
      std::sample(samples.begin(), std::next(samples.begin(),fittest), std::back_inserter(parents), 2, gen);
      *it = Scheme::generate(parents[0],parents[1]);
      it->mutate(0.15, fitness);
    }
  }
}
