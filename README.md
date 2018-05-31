
# <img src="https://raw.githubusercontent.com/ITBE-Lab/ma/release/MA.png" align="center" width="90"> The Modular Aligner
MA is a novel open source tool for the efficient and accurate alignment of short and long reads of various sequencers. The approach has a highly modular architecture and everyone is invited to propose/integrate 
new modules for getting better results overall or for specific sequencers. The
design aims at a smooth Python integration, while keeping the performance delivered by C++. So, the
general idea is the coupling of several C++-modules via Python, where each module performs a different 
part of the alignment  process. Well performing module combinations are finally coupled
under the roof of a single C++ application that embodies a single aligner 
for high throughput computing.
## Algorithmic Approaches
Many long and short read aligners conceptually follow a 3 stage approach:
1. **Seeding**. Using some form of indexing mechanism a set of perfect matches (or with few mismatches) is computed for the given query.
2. **Coupling**. Identification of promising seed subsets that could lead to an optimal alignment. A 
popular technique for this task is *chaining*.
3. **Dynamic Programming**. Filling of remaining gaps between seeds and extension of the outer 
end points given by the seeds.

MA stays with this basic pattern but delivers novel algorithmic solutions for the stages 1 and 2. 
MA introduces a divide and conquer approach for seeding on the Foundation of the FMD-Index. The 
advantage of this seeding variant is the reduction of the overall number of seeds compared to a
classical FMD-index based extension process (e.g. in BWA-MEM). One could expect 
that a reduction of the number seeds hurts an aligner's accuracy, but this 
is only true, if we loose significant and insignificant seeds at equal rates. And, as this latter 
statement already suggests, the divide and conquer variant goes after the significant seeds.\
In stage 2 we come up with two new ideas. First, the *strip of consideration* (SoC) and, 
second, *seed harmonization*. Both techniques represent simple line-sweeps that iterate over all seeds 
(or some already heuristically filtered set of seeds). The SoC is used to quickly identify 
promising regions on the reference, while the seed harmonization purges contradicting 
seeds by relying on a virtual guide line. The highlights of these approaches are: No specially
tailored data structures required. Highly efficient. Easy to implement (roughly 50 lines Pseudocode).

## Python Integration

> **Python Integration is an optional feature of MA.** By default MA is built without Python support.

The Python integration of MA is done via [Boost.Python](https://www.boost.org/ "Boost.Python"). 
So, it is necessary to have a boost deployment with Python3 support. If there is a Boost.Python 
library, then the Makefile can be used with `WITH_PYTHON=1`. The folder *MA* of the repository 
contains a Python3 module that incorporates MA as Python3 module. Please note that the Python3 module relies on
the shared library *libMA.so* that is built by the Makefile. Alignments using the Python integration 
can be done without a computational penalty. Python is only responsible for an initial C++-module 
coupling, while all actual computations are done within the C++-modules. The idea is similar 
to the one used in the context of [TensorFlow](https://www.tensorflow.org "TensorFlow").

## Getting Started
MA is currently still in an evaluation and testing process. Feedback about bugs is highly welcomed. 
The below procedure was checked on a fresh *Debian 4.9.88* with the following additional 
packages installed: `git`, `make` and `build-essential`.

 
 Get the github clone and call make. (Available make switches are documented below):

    git clone https://github.com/ITBE-Lab/ma
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

## Makefile switches


- WITH_AVX:
    Enables avx instructions for unused part of code (the sw implementation)
    May cause invalid instruction errors if the hardware does not support AVX2

- WITH_PYTHON:
    Compiles MA as a shared library which can be imported by python.
    The ma executable is then dynamically linked to the shared library, 
    so the current folder must be added to LD_LIBRARY_PATH.
    This options requires that $(BOOST_ROOT), $(PYTHON_INCLUDE) and $(PYTHON_LIB) are correctly set.
    *Tested using boost 1.65.1 (with the boost python library) and python 3.*

- DEBUG:
    Compiles the code un-optimized with assertions enabled and multiple self-checks during runtime. 
    Mostly intended for debugging purposes.

## MA options

```
General options:
    -h, --help              Display the complete help screen
    -a, --align             Do sequence alignment
    -t, --threads arg       Used concurrency
    -f, --fmdIndex          Do FMD-index generation

Alignment options (requires -a):
    -i, --alignIn args      Input file(s) as (multi-)fasta(-q)
    -o, --alignOut arg      Output file as SAM
    -g, --genome arg        FMD-index input file prefix
    -p, --parameterSet arg  Pre-setting [fast/accurate]
    -s, --seedSet arg       Seeding strategy [SMEMs/maxSpanning]
    -n, --reportN arg       Report <= N alignments; 0: unlimited
    -l, --minLen arg        Minimum seed length
    -v, --giveUp arg        Minimum SoC score (relative to query length)
    -b, --basicMode         Disable DP
        --Match arg         DP match score.
        --MissMatch arg     DP mismatch penalty.
        --Gap arg           DP gap open penalty.
        --Extend arg        DP gap extend penalty.

Paired Reads options (requires either -U or -N):
    -U, --uniform           Enable paired alignment; Distance as uniform distribution
    -N, --normal            Enable paired alignment; Distance as normal distribution
    -u, --unpaired arg      Penalty for unpaired alignments
    -m, --mean arg          Gap distance mean
    -d, --std arg           Gap distance standard deviation

FMD-Index Generation options (requires -f):
    -I, --indexIn args      (Multi-)Fasta input file(s)
    -O, --indexOut arg      FMD-index output file prefix
```


## Thanks ...

MA relies on the hard work of other projects. These are:

- Dynamic programming with adaptive band: https://github.com/ocxtal/libgaba
- Dynamic programming with static band: https://github.com/lh3/ksw2
- Cmd line option parser: https://github.com/jarro2783/cxxopts
- Parts of the code for FMD-index creation and extension were picked from https://github.com/lh3/bwa


Many thanks to the creators and contributors of the above projects ...

## Authors

MA is created and maintained by Markus Schmidt and Arne Kutzner.
