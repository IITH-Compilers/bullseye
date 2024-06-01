# Bullseye

The BullsEye cache model is an analytical tool that calculates the number of cache misses for a given affine program.
BullsEye is developed by the IITH Compilers group, and the cache model is based on implementation from HayStack (https://dl.acm.org/doi/10.1145/3314221.3314606 appeared in PLDI 2019).

## Installation

### Requirements

Before installing the package make sure you install the following dependencies:
- llvm (version 10)
- clang (version 10)
- gmp
- ntl
- boost (program options)
- libyaml
- cplex (Installed in $CPLEX_HOME(Example: /opt/CPLEX_Studio201/))

If you are using an Ubuntu 20.04 system, you can install these dependencies by executing the following command:

```bash
sudo apt-get install llvm-10-dev libclang-10-dev libgmp3-dev libntl-dev libboost-program-options-dev libyaml-dev
```
### Build the Project

#### Clone the Repository
The following command will clone this repository and initialize all the submodules.
```bash
git clone git@github.com:IITH-Compilers/bullseye.git
cd bullseye 
git submodule update --init --recursive
```
#### Setting thing Up
```bash
./isl_fix.sh
./get_submodules.sh
./autogen.sh
```
These commands will apply the isl fix and initialize Autotools. With these steps completed, you can now configure and build.

#### Configure build
```bash
export CPLEX_HOME=path_to_CPLEX
./configure --prefix=$HOME LIBS='-ldl'
make
make install
```
First export the path of CPLEX installation and then build the project.
By setting the prefix to the home folder, we can choose to install the tool in the user directory, although this is optional. Additionally, we can specify and configure options to direct Autotools to dependencies like the boost program options library if they are not detected automatically.

> [!WARNING]
> If you have multiple versions of llvm and or clang installed you might have multiple
> llvm-config, in such a case manually give a path to the llvm-config-10 as follows
> ```bash
> export LLVM_CONFIG=path_to_llvm-config-10
> export CPLEX_HOME=path_to_CPLEX
> ./configure --prefix=$HOME LIBS='-ldl'
> make
> make install 
> ```

## Usage

Use the following command to run bullseye on the example program gemm.c
```
bullseye -f ./examples/gemm.c
```

## Citing BullsEye

If you use BullsEye, please consider citing the relevant publications:

``` bibtex
@article{10.1145/3558003,
author = {Shah, Nilesh Rajendra and Misra, Ashitabh and Min\'{e}, Antoine and Venkat, Rakesh and Upadrasta, Ramakrishna},
title = {
BullsEye : Scalable and Accurate Approximation Framework for Cache Miss Calculation},
year = {2022},
issue_date = {March 2023},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {20},
number = {1},
issn = {1544-3566},
url = {https://doi.org/10.1145/3558003},
doi = {10.1145/3558003},
abstract = {For Affine Control Programs or Static Control Programs (SCoP), symbolic counting of reuse distances could induce polynomials for each reuse pair. These polynomials along with cache capacity constraints lead to non-affine (semi-algebraic) sets; and counting these sets is considered to be a hard problem. The state-of-the-art methods use various exact enumeration techniques relying on existing cardinality algorithms that can efficiently count affine sets.We propose BullsEye , a novel, scalable, accurate, and problem-size independent approximation framework. It is an analytical cache model for fully associative caches with LRU replacement policy focusing on sampling and linearization of non-affine stack distance polynomials. First, we propose a simple domain sampling method that can improve the scalability of exact enumeration. Second, we propose linearization techniques relying on Handelman’s theorem and Bernstein’s representation. To improve the scalability of the Handelman’s theorem linearization technique, we propose template (Interval or Octagon) sub-polyhedral approximations.Our methods obtain significant compile-time improvements with high-accuracy when compared to HayStack on important polyhedral compilation kernels such as nussinov, cholesky, and adi from PolyBench, and harris, gaussianblur from LLVM-TestSuite. Overall, on PolyBench kernels, our methods show up to 3.31\texttimes{} (geomean) speedup with errors below ≈ 0.08\% (geomean) for the octagon sub-polyhedral approximation.},
journal = {ACM Trans. Archit. Code Optim.},
month = {nov},
articleno = {2},
numpages = {28},
keywords = {Static analysis, cache model, performance analysis}
}
```
