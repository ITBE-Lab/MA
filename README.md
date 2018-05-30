
## Getting Started

execute following commands:

    git clone https://github.com/ItBeLab/ma
    cd ma
    make

Test your installation with:

    ./ma -h

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


