#include "genetic.h"
#include <list>
#include <random>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <execution>
// #include <pstl/execution>
// #include <pstl/algorithm>

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
      [=](double val){return lhs*val;});
  return result;
}

struct ScalarTimesVector {
  VecD operator()(const double& lhs, const VecD& rhs) {
    return lhs*rhs;
  }
};

struct VectorPlusVector {
  VecD operator()(VecD const& lhs, VecD const& rhs) {
    return lhs+rhs;
  }
};

Scheme::Scheme(int n, double dt) {
  this->n = n;

  a_lower = VecD(static_cast<int>(n*(n-1)/2));
  b = VecD(n);
  c = VecD(n);
  // a = std::vector<VecD>(n,VecD(n));
  this->dt = dt;
}

void Scheme::init() {
  std::random_device rd{};
  std::mt19937 gen(rd());
  std::normal_distribution<> gaussian(0.5,0.5); // (Mean, Standard deviation)

  auto gaussian_gen = std::bind(gaussian,gen);
  std::generate(a_lower.begin(),a_lower.end(),gaussian_gen);
  std::generate(b.begin(),b.end(),gaussian_gen);
  std::generate(c.begin(),c.end(),gaussian_gen);
  /* dt = std::pow(10,-3*gaussian_gen()); */
}

void Scheme::mutate(double prob, double scaling) {
  std::random_device rd{};
  std::mt19937 gen(rd());
  std::normal_distribution<> gaussian(0,3); // (Mean, Standard deviation)
  std::normal_distribution<> scaled_gaussian(0.,std::pow(scaling,2)); // (Mean, Standard deviation)

  // Only with probability prob
  std::bernoulli_distribution decider(prob);
  std::bernoulli_distribution decider5050(1./2);
  while (decider(gen)) {
    // select one property
    std::vector<VecD*> l = {&this->a_lower, &this->b, &this->c};
    std::vector<VecD*> chosen_property;
    std::sample(l.begin(), l.end(), std::back_inserter(chosen_property),
        1, gen);
    if (decider5050(gen))
      chosen_property[0]->at(std::rand()%chosen_property[0]->size()) += scaled_gaussian(gen);
    else
      chosen_property[0]->at(std::rand()%chosen_property[0]->size()) = gaussian(gen);
  }
}

Scheme Scheme::generate(const Scheme& lhs, const Scheme& rhs) {
  std::random_device rd{};
  std::mt19937 gen(rd());
  std::bernoulli_distribution decider_one_third(1./3);
  std::bernoulli_distribution decider_one_half(1./2);

  Scheme s(lhs.n, lhs.dt);
  if (lhs.n != rhs.n) throw "Uneqal sizes n!";
  s.a_lower = {};
  s.b = {};
  s.c = {};
  std::transform(lhs.a_lower.begin(), lhs.a_lower.end(), rhs.a_lower.begin(),
      std::back_inserter(s.a_lower), [&](double v1, double v2) {if (decider_one_third(gen)) return v1; else if (decider_one_half(gen)) return v2; else return (v1+v2)/2.;});
  std::transform(lhs.b.begin(), lhs.b.end(), rhs.b.begin(),
      std::back_inserter(s.b),       [&](double v1, double v2) {if (decider_one_third(gen)) return v1; else if (decider_one_half(gen)) return v2; else return (v1+v2)/2.;});
  std::transform(lhs.c.begin(), lhs.c.end(), rhs.c.begin(),
      std::back_inserter(s.c),       [&](double v1, double v2) {if (decider_one_third(gen)) return v1; else if (decider_one_half(gen)) return v2; else return (v1+v2)/2.;});
  s.dt = (lhs.dt+rhs.dt)/2.;
  return s;
}

std::vector<VecD> Scheme::run(const VecD& x0, std::function<VecD(double,VecD)> f, double t_end) const {
  int n_t = std::floor((t_end-0)/this->dt)+1;
  std::vector<VecD> x(n_t,VecD(x0.size(),0));
  x[0] = x0;
  double t = 0; // Starting at fist timestep
  std::vector<VecD> k(this->n,VecD(x0.size(),0.));
  for (int i = 0; i < n_t-1; i++) {
    k = this->calcKVec(t,x[i],f);
    x[i+1] = std::transform_reduce(
      std::execution::unseq,
      this->b.begin(), this->b.end(),
      k.begin(), x[i], 
      VectorPlusVector(), ScalarTimesVector()
    );
    /*
    std::vector<VecD> b_k; b_k.reserve(k.size());
    std::transform(b.begin(),b.end(),
                   k.begin(),std::back_inserter(b_k),
                   [](double b, const VecD& k) {return b*k;});
    x[i+1] = x[i] + this->dt * std::accumulate(b_k.begin(),b_k.end(),VecD(x[i].size(),0.));
    t += this->dt;
    */
    /* x += this->dt * std::transform_reduce(this->b.begin(), this->b.end(), k.begin(), VecD(x0.size(),0.), */
    /*                                       [](double b, VecD k){return b*k;},[](const VecD& l, const VecD& r){return l+r;}); */
  }
  return x;
}

std::vector<VecD> A_full_from_lower(VecD const& a_lower) {
  size_t dims = static_cast<size_t>(std::sqrt(0.25+2*a_lower.size())-1./2) + 1;  // a_lower.size() = n(n+1)/2 => pq-equation: n = -1/2 + sqrt(1/4 + a_lower.size)
  std::vector<VecD> A(dims, VecD());
  size_t lower_idx = 0;

  for (size_t row = 1; row < dims; row++) {
    A[row] = VecD(row);
    for (size_t col = 0; col < row; col++) {
      A[row][col] = a_lower[lower_idx];
      lower_idx++;
    }
  }

  return A;
}

std::vector<VecD> Scheme::calcKVec(double t, const VecD& x, std::function<VecD(double,VecD)> f) const {
  std::vector<VecD> k(this->n,VecD(x.size(),0.));
  std::vector<VecD> A = A_full_from_lower(this->a_lower);
  double t_inner;
  VecD x_inner;
  for (size_t j = 0; j < this->n; j++) {
    t_inner = t + this->dt*this->c[j];
    x_inner = x;
    for (size_t l = 0; l < j; l++) {
      x_inner += this->dt * A[j][l]*k[l];
    }
    k[j] = f(t_inner, x_inner);
  }

  return k;
}

/**
std::vector<VecD> Scheme::calcKVec(double t, const VecD& x, std::function<VecD(double,VecD)> f) const {
  int a_idx, k_idx, jump;
  double t_inner;
  std::vector<VecD> k(this->n,VecD(x.size(),0.));
  VecD x_inner;
  for (int j = 0; j < this->n; j++) {
    t_inner = t + this->c[j]*this->dt;
    x_inner = x;
    k_idx = 0;
    a_idx = j-1;
    jump = this->n-2;
    while (k_idx < j) {
      x_inner += this->dt*this->a_lower[a_idx]*k[k_idx];
      a_idx += jump;
      jump--;
      k_idx++;
    }
    k[j] = f(t+this->dt*this->c[j], x_inner);
  }
  return k;
}
*/