#include <utility>

#pragma once

#define CHECK_VAR_INT_RANGE(VAR, LOW, HIGH) if (!this->check_int_range(VAR, LOW, HIGH, #VAR)) return;
#define CHECK_VAR_REAL_RANGE(VAR, LOW, HIGH) if (!this->check_real_range(VAR, LOW, HIGH, #VAR)) return;
#define CHECK_VAR_REAL_RANGE_2(VAR, LOW, HIGH) if (!this->check_real_range(VAR, LOW, HIGH, #VAR)) return VAR;

template<typename Real, typename Complex>
class Integral {
private:
    void print_error(const char *name, const Real &val) {
        fmt::println(std::cerr,
                     "[{}] a = {}, r_s = {}, theta_s {}, r_o = {}, lambda = {}, eta = {}, {} = {}, out of range",
                     child_class_name,
                     data.a,
                     data.r_s,
                     data.theta_s,
                     data.r_o,
                     data.lambda,
                     data.eta, name, val);
    }

public:
    ForwardRayTracing<Real, Complex> &data;
    std::string child_class_name;

    explicit Integral(ForwardRayTracing<Real, Complex> &data_, std::string child_class_name_) : data(data_),
                                                                                                child_class_name(
                                                                                                        std::move(
                                                                                                                child_class_name_)) {}

    bool check_int_range(const Real &val, int low, int high, const char *name) {
        if (val < low || val > high) {
            print_error(name, val);
            data.ray_status = RayStatus::INTERNAL_ERROR;
            return false;
        } else {
            return true;
        }
    }

    bool check_real_range(const Real &val, const Real &low, const Real &high, const char *name) {
        if (val < low || val > high) {
            print_error(name, val);
            data.ray_status = RayStatus::INTERNAL_ERROR;
            return false;
        } else {
            return true;
        }
    }
};