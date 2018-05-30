
## Getting Started

execute following commands in order to **install MA**:

    git clone https://github.com/ItBeLab/ma
    git checkout v0.1.0-alpha # (*)
    cd ma
    make

(*) This command is used in order to switch to the last stable release. 
Skip it if you want to use the newest available version.

Then test your installation with:

    ./ma --help

This should display the help screen of MA.

\
In order to create alignments you need to **compute an index** for your reference genome first.
This is done using following command:

    ./ma --fmdIndex --indexIn <filename_of_genome> --indexOut <filename_of_index>

Where `<filename_of_genome>` if the path to a (multi-)fasta file containing the reference genome and 
`<filename_of_index>` is the filename-prefix that shall be used to store the index.

Now you can align using following command in order to **align**:

    ./ma --align --genome <filename_of_index> --alignIn <fasta_in> --alignOut <sam_out>

Where `<fasta_in>` is a (multi-)fasta(-q) file containing the queries.
`<sam_out>` is the filename of the output file that shall be created.
You can switch between **MA accurate** and **MA fast** by using the `--parameterset accurate` or 
`--parameterset fast` switch, respectively.

---
Tested using a vanilla Debian 4.9.88 with following additional packages installed:
`git`, `make` and `build-essential`.

## Citing MA

MA is unpublished so far.

## Compilation switches

[\
    WITH_AVX:
        Enables avx instructions for unused part of code (the sw implementation)
        MAy cause invalid instruction errors if the hardware does not support AVX2\
]

- WITH_PYTHON:
    Compiles MA as a shared library which can be imported by python.
    The ma executable is then dynamically linked to the shared library, 
    so the current folder must be added to LD_LIBRARY_PATH.
    This options requires that $(BOOST_ROOT), $(PYTHON_INCLUDE) and $(PYTHON_LIB) are correctly set.
    *Tested using boost 1.65.1 (with the boost python library) and python 3.*

- WITH_GPU_SW:
    Compiles the gpu sw implementation.
    In combination with WITH_PYTHON the gpu_sw is also accessible via python.
    *Requires a cuda compiler.* 

- DEBUG:\
    Compiles the code un-optimized with assertions enabled and multiple self-checks during runtime. 
    Mostly intended for debugging purposes.

## Thanks to

We integrated several other projects (some only in parts).
Here are their github pages:

- Dynamic programming with adaptive band: https://github.com/ocxtal/libgaba
- Dynamic programming with static band: https://github.com/lh3/ksw2
- Cmd line option parser: https://github.com/jarro2783/cxxopts
- Initially we started on the basis of bwa-mem; 
    However, by now, we replaced almost everything with our own code. 
    Still, the FMD-index, Pack and SMEM-extension remain highly similar to those found here: https://github.com/lh3/bwa
- Ransac implementation ?

Many thanks to the creators of above projects. 
And special thanks to Li Heng, whose aligner served as our starting point.
