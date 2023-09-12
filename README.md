# ForwardRayTracing
 
```yaml
spack:
  specs:
    - boost+filesystem+program_options
    - catch2
    - fmt
    - gmp
    - mpfr
    - mpc
    - eigen
    - intel-oneapi-tbb
    - python
    - py-pybind11
  view: true
  concretizer:
    unify: true
```

```vcpkg
git clone https://github.com/Microsoft/vcpkg.git
.\vcpkg\bootstrap-vcpkg.bat
.\vcpkg\vcpkg.exe install boost catch2 fmt eigen3 tbb pybind11
```