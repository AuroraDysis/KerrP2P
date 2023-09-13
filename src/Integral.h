#include <utility>

#pragma once

#define CHECK_STATUS if (this->data.ray_status != RayStatus::NORMAL) return;
#define CHECK_VAR(VAR, COND) if (!this->check_variable(VAR, COND, #VAR)) return;

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

    bool check_variable(const Real &val, bool condition, const char *name) {
        if (!condition) {
            print_error(name, val);
            data.ray_status = RayStatus::INTERNAL_ERROR;
            return false;
        } else {
            return true;
        }
    }
};