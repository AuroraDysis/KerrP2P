#pragma once

#define CHECK_SIN_RANGE(VAR) if (!this->check_sin_range(VAR, #VAR)) return;
#define CHECK_COS_RANGE(VAR) if (!this->check_cos_range(VAR, #VAR)) return;

template <typename Real, typename Complex>
class Integral {
public:
	ForwardRayTracing<Real, Complex>& data;

	explicit Integral(ForwardRayTracing<Real, Complex>& data_) : data(data_) {}

	bool check_sin_range(const Real &sin_x, const char *name) {
		if (sin_x < -1 || sin_x > 1) {
			fmt::println(std::cerr, "lambda = {}, eta = {}, {} = {}, out of range", data.lambda, data.eta, name, sin_x);
			data.ray_status = RayStatus::INTERNAL_ERROR;
			return false;
		} else
		{
			return true;
		}
	}

	bool check_cos_range(const Real& cos_x, const char* name) {
		if (cos_x < -1 || cos_x > 1) {
			fmt::println(std::cerr, "lambda = {}, eta = {}, {} = {}, out of range", data.lambda, data.eta, name, cos_x);
			data.ray_status = RayStatus::INTERNAL_ERROR;
			return false;
		}
		else
		{
			return true;
		}
	}
};