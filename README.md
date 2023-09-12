# ForwardRayTracing

## Quick Start

### Linux

[Spack](https://spack.io/)

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

### MacOS

[Homebrew](https://brew.sh/)

```bash
brew install boost catch2 fmt gmp mpfr mpc eigen tbb python pybind11
```

### Windows ()

[vcpkg](https://vcpkg.io/en/index.html)

```bash
vcpkg.exe install boost catch2 fmt eigen3 tbb python3 pybind11
```

## Acknowledgements

