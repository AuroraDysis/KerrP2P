#include <utility>

#pragma once

#define CHECK_VAR_INT_RANGE(VAR, LOW, HIGH) if (!this->check_int_range(VAR, LOW, HIGH, #VAR)) return;

template<typename Real, typename Complex>
class Integral {
private:
    void print_error(const char *name, const Real &x) {
        fmt::println(std::cerr,
                     "[{}] a = {}, r_s = {}, theta_s {}, r_o = {}, lambda = {}, eta = {}, {} = {}, out of range",
                     child_class_name,
                     data.a,
                     data.r_s,
                     data.theta_s,
                     data.r_o,
                     data.lambda,
                     data.eta, name, x);
    }

public:
    ForwardRayTracing<Real, Complex> &data;
    std::string child_class_name;

    explicit Integral(ForwardRayTracing<Real, Complex> &data_, std::string child_class_name_) : data(data_),
                                                                                                child_class_name(
                                                                                                        std::move(
                                                                                                                child_class_name_)) {}

    bool check_int_range(const Real &x, int low, int high, const char *name) {
        if (x < low || x > high) {
            print_error(name);
            data.ray_status = RayStatus::INTERNAL_ERROR;
            return false;
        } else {
            return true;
        }
    }

    bool check_real_range(const Real &x, const Real &low, const Real &high, const char *name) {
        if (x < low || x > high) {
            print_error(name);
            data.ray_status = RayStatus::INTERNAL_ERROR;
            return false;
        } else {
            return true;
        }
    }
};