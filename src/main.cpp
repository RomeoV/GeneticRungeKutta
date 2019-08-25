#include <cmath>
#include <iostream>
#include <algorithm>
#include <random>
#include <numeric>
#include <iomanip>
#include <execution>
#include <map>
#include "genetic.h"

using VecS = std::vector<Scheme>;

std::pair<double,VecD> evalSinCos(const Scheme& s) {
  double t_end = M_PI;
  auto [timesteps, result] = s.run({0.,1.},[](double t, const VecD& v)->VecD {return {v[1],-v[0]};}, t_end);
  double error = 0;
  size_t n_t = timesteps.size();
  double t = 0;
  for (size_t i = 0; i < n_t; i++) {
    error += norm_dist<2>(result[i], {std::sin(t), std::cos(t)});
    t += s.dt;
  }

  if (std::isnan(error)) {
    std::cerr << "Evaluation result is NaN!" << std::endl;
    std::cerr << s << std::endl;
    throw std::runtime_error("Evaluation result is NaN!");
  }
  return {error, result.back()};
}

std::pair<double,VecD> eval2DSecondOrder(const Scheme& s) {
  double t_end = 0.1;
  auto [timesteps, result] = s.run({1.,1.},[](double t, const VecD& v)->VecD {return {v[1],-10*v[1]+std::pow(t,2)};}, t_end);
  double error = 0;

  size_t n_t = timesteps.size();
  double t = 0;
  for (size_t i = 0; i < n_t; i++) {
    double C1 = -499./5000;
    double C2 = 1-C1;
    double y = std::pow(t,3)/30. - std::pow(t,2)/.100 + t/500. + C1*std::exp(-10.*t) + C2;
    error += std::pow(result[i][0]-y,2);
    t += s.dt;
  }

  return {error, result.back()};
}


int main() {
  std::random_device rd{};
  std::mt19937 gen(rd());
  const int sample_size = 100;
  const int epochs = 500;
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
    std::map<Scheme, double> fitness_map;
    std::transform(
        samples.begin(), samples.end(),
        std::inserter(fitness_map, fitness_map.begin()),
        [](Scheme const& s)->std::pair<Scheme, double> {
          try {return {s, evalSinCos(s).first};}
          catch (std::runtime_error const& e) {
              return {s, std::numeric_limits<double>::max()};
          }}
     );
    std::sort(//std::execution::par,
        samples.begin(), samples.end(), 
        [&fitness_map](const Scheme& s1, const Scheme& s2){
          return fitness_map[s1] < fitness_map[s2];
        });
    std::clog << "Done." << std::endl;
    std::clog << "Evaluating best sample..." << std::endl;
    auto [best_fitness, x_end] = evalSinCos(samples[0]);
    std::clog << std::setprecision(15);
    std::clog << "Best Score was: " << best_fitness << "!" << std::endl;
    std::clog << "x = " << x_end << std::endl;
    std::clog << std::setprecision(5);
    std::clog << samples[0] << std::endl;
    std::clog << "#####################" << std::endl << std::endl;;

    // replace unfit by offspring of fittest population
    for (auto it = std::next(samples.begin(),fittest); it != samples.end(); it++) {
      std::vector<Scheme> parents;
      std::sample(samples.begin(), std::next(samples.begin(),fittest), std::back_inserter(parents), 2, gen);
      *it = Scheme::generate(parents[0],parents[1]);
      it->mutate(0.15, best_fitness);
    }
  }
}
