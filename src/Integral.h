#pragma once

#define CHECK_VAR_INT_RANGE(VAR, LOW, HIGH) if (!this->check_int_range(VAR, LOW, HIGH, #VAR)) return;

template <typename Real, typename Complex>
class Integral {
public:
	ForwardRayTracing<Real, Complex>& data;
	std::string child_class_name;

	explicit Integral(ForwardRayTracing<Real, Complex>& data_, std::string child_class_name_) : data(data_), child_class_name(child_class_name_) {}

	bool check_int_range(const Real& x, int low, int high, const char* name) {
		if (x < low || x > high) {
			fmt::println(std::cerr, "[{}] lambda = {}, eta = {}, {} = {}, out of range", child_class_name, data.lambda, data.eta, name, x);
			data.ray_status = RayStatus::INTERNAL_ERROR;
			return false;
		}
		else
		{
			return true;
		}
	}

	bool check_real_range(const Real& x, const Real& low, const Real& high, const char* name) {
		if (x < low || x > high) {
			fmt::println(std::cerr, "[{}] lambda = {}, eta = {}, {} = {}, out of range", child_class_name, data.lambda, data.eta, name, x);
			data.ray_status = RayStatus::INTERNAL_ERROR;
			return false;
		}
		else
		{
			return true;
		}
	}
};