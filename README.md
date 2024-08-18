# KerrP2P

`KerrP2P` is a software package designed for forward ray tracing in Kerr spacetime. It is specifically tailored to calculate **multiple** null geodesics between designated "source" and "observer" points as detailed in the paper 'Forward Ray Tracing and Hot Spots in Kerr Spacetime' by Lihang Zhou, Zhen Zhong, Yifan Chen, and Vitor Cardoso (The preprint will be available soon on arXiv).

The software provides precise solutions to the geodesic equations using Jacobi elliptic functions, based on the methodology outlined in the paper 'The Null Geodesics of the Kerr Exterior' ([arXiv:1910.12881](https://arxiv.org/abs/1910.12881)).

[![Paper](https://img.shields.io/badge/Paper-arXiv-blue)](https://arxiv.org/abs/your_paper_link)
[![License](https://img.shields.io/badge/License-MIT-green)](https://opensource.org/licenses/MIT)

# Features

- Forward ray tracing in Kerr spacetime: this involves calculating **multiple** null geodesics from a given source to an observer
- Support for arbitrary precision arithmetic
- Python bindings for easy to use interface
- Support for multiple platforms: Linux, macOS, and Windows

## Installation

### Prerequisites

Before installing `KerrP2P`, ensure that you have the following dependencies installed on your system:

- Boost (with filesystem and program_options components)
- Catch2 (optional, for testing)
- fmt
- GMP (optional)
- MPFR (optional)
- MPC (optional)
- Eigen
- Intel oneAPI TBB
- Python (optional)
- pybind11 (optional)

### Linux

On Linux, you can use [Spack](https://spack.io/) to install the dependencies:

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

### macOS

On macOS, you can use [Homebrew](https://brew.sh/) to install the dependencies:

```bash
brew install boost catch2 fmt gmp mpfr mpc eigen tbb python pybind11
```

### Windows

On Windows, you can use [vcpkg](https://vcpkg.io/en/index.html) to install the dependencies:

```bash
vcpkg.exe install boost catch2 fmt eigen3 tbb python3 pybind11
```

## Usage

To use `KerrP2P`, follow these steps:

1. Clone the repository:

```bash
git clone https://github.com/your_username/KerrP2P.git
```

2. Build the project using CMake.

```bash
cd KerrP2P
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Contributing

Contributions to `KerrP2P` are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request on the [GitHub repository](https://github.com/AuroraDysis/KerrP2P).

## License

`KerrP2P` is released under the [MIT License](https://opensource.org/licenses/MIT). See the [LICENSE](LICENSE) file for more details.

## Citation

If you use `KerrP2P` in your research, please cite the following paper:

> [!CAUTION]
> The preprint will be available soon on arXiv.

```bibtex
@article{KerrP2P,
  title = "{Forward Ray Tracing and Hot Spots in Kerr Spacetime}",
  author = "Lihang Zhou, Zhen Zhong, Yifan Chen and Vitor Cardoso",
  year = "2024",
}
```
