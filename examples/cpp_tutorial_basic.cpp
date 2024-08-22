#include "ForwardRayTracing.h"
#include "Utils.h"

#include <fmt/ranges.h>

using std::string;

int main(int argc, char *argv[]) {
  using Real = Float256;
  using Complex = Complex256;
  ForwardRayTracingParams<Real> params;

  const auto &pi = boost::math::constants::pi<Real>();
  params.a = boost::lexical_cast<Real>("0.8");
  params.r_s = boost::lexical_cast<Real>("1.7");
  params.theta_s = boost::lexical_cast<Real>("89.9") * pi / 180;
  params.r_o = 1000;
  params.nu_r = Sign::POSITIVE;
  params.nu_theta = Sign::NEGATIVE;

  // Set parameters from rc, log_abs_d, d_sign
  // params.rc = boost::lexical_cast<Real>("1.8117208808167675");
  // params.log_abs_d = boost::lexical_cast<Real>("0.34917458729364625");
  // params.d_sign = Sign::NEGATIVE;
  // params.rc_d_to_lambda_q();

  // set parameters from lambda, q
  params.lambda = boost::lexical_cast<Real>("-0.3991390432493585887150804403");
  params.q = sqrt(boost::lexical_cast<Real>("4.886104404637876423649492482"));
  params.calc_t_f = true;

  auto forward = ForwardRayTracing<Real, Complex>::get_from_cache();
  forward->calc_ray(params);
  fmt::println("ray status: {}", ray_status_to_str(forward->ray_status));
  fmt::println("theta_f: {}, phi_f: {}", forward->theta_f, forward->phi_f);
}
