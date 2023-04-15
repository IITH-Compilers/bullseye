# Bullseye

The BullsEye cache model is an analytical tool that calculates the number of cache misses for a given affine program.
BullsEye is developed by IITH Compilers group and adopted from HayStack(PLDI 2019).

## Installation

Before installing the package make sure you install die following dependencies:
- llvm/clang
- gmp
- ntl
- boost (program options)
- libyaml
- cplex (Installed in $CPLEX_HOME(Example: /opt/CPLEX_Studio201/))

If you are using an Ubuntu 20.04 system, you can install these dependencies by executing the following command:
```
sudo apt-get install llvm-dev libclang-dev libgmp3-dev libntl-dev libboost-program-options-dev libyaml-dev
```
After installing all the dependencies, navigate to the Haystack folder and execute the following commands:
```
./get_submodules.sh
./autogen.sh
```
These commands will fetch all the submodules and initialize autotools. With these steps completed, you can now configure and build.
First export the path of CPLEX installation and then build the project.
```
export CPLEX_HOME=path_to_CPLEX
./configure --prefix=$HOME
make
make install
```
By setting the prefix to the home folder, we can choose to install the tool in the user directory, although this is optional. Additionally, we can specify configure options to direct autotools to dependencies like the boost program options library, in case they are not detected automatically.

## Usage

Use the following command to run bullseye on example program gemm.c
```
bullseye -f ./examples/gemm.c
```  

