
## Getting Started

execute following commands:

    git clone https://github.com/ItBeLab/ma
    git checkout v0.1.0-alpha # (*)
    cd ma
    make

(*) Execute this command in order to go to the last stable release
Test your installation with:

    ./ma -h

You should see the help menu of ma.

Tested using a vanilla Debian 4.9.88 with following additional packages installed:
git, make and build-essential.

## Citing MA

MA is unpublished so far.

## Compilation switches

[
    WITH_AVX
        Enables avx instructions for unused part of code (the sw implementation)
        MAy cause invalid instruction errors if the hardware does not support AVX2
]

WITH_GPU_SW
    Compiles the gpu sw implementation.
    In combination with WITH_PYTHON the gpu_sw is also accessible via python

WITH_PYTHON
    Compiles MA as a shared library libMA.so which can be imported by python.
    The ma executable is then dynamically linked to the shared library, 
    so the current folder must be in the system variable.
    This requires that $(BOOST_ROOT), $(PYTHON_INCLUDE) and $(PYTHON_LIB) are correctly set.
    Tested using boost 1.65.1 (with the boost python library) and python 3.

DEBUG
    Compiles the code un-optimized with assertions enabled and multiple runtime self-checks.

## Thanks to

We integrated several other projects (some only in parts).
Here are their github pages:

- https://github.com/ocxtal/libgaba (dynamic programming with adaptive band)
- https://github.com/lh3/ksw2 (dynamic programming with static band)
- https://github.com/jarro2783/cxxopts (cmd line option parser)
- https://github.com/lh3/bwa (initially we started on the basis of bwa-mem; 
    however by now we replaced almost everything with our own code. 
    Still, the FMD-index, Pack and the SMEM-extension remain highly similar to those found here.)
- ransac implementation ?

