# Bullseye

The BullsEye cache model is an analytical tool that calculates the number of cache misses for a given affine program.
BullsEye is developed by IITH Compilers group and adopted from HayStack(PLDI 2019).

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
The following command will clone this repository initialize all the submodules.
```bash
git clone git@github.com:IITH-Compilers/bullseye.git
cd byllseye
git checkout libbullseye
git pull
git submodule update --init --recursive
```
#### Setting thing Up
```bash
./isl_fix.sh
./autogen.sh
```
These commands will apply the isl fix and initialize autotools. With these steps completed, you can now configure and build.

#### Configure build
```bash
export CPLEX_HOME=path_to_CPLEX
./configure --prefix=$HOME LIBS='-ldl'
make
make install
```
First export the path of CPLEX installation and then build the project.
By setting the prefix to the home folder, we can choose to install the tool in the user directory, although this is optional. Additionally, we can specify configure options to direct autotools to dependencies like the boost program options library, in case they are not detected automatically.

> [!WARNING]
> If you have multiple versions of llvm and or clang installed you might have muliple
> llvm-config, in such a case manually give path to the llvm-config-10 like follows
> ```bash
> export LLVM_CONFIG=path_to_llvm-config-10
> export CPLEX_HOME=path_to_CPLEX
> ./configure --prefix=$HOME LIBS='-ldl'
> make
> make install 
> ```

## Usage

Use the following command to run bullseye on example program gemm.c
```
bullseye -f ./examples/gemm.c
```  

