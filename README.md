
# <img src="https://raw.githubusercontent.com/ITBE-Lab/MA/release/MA.png" align="center" width="90"> The Modular Aligner
MA is a novel open source application for the efficient and accurate alignment of short and long reads of various sequencers. The approach has a highly modular architecture and everyone is invited to propose/integrate new modules. 
The design aims at a smooth Python integration, while keeping the performance delivered by C++. So, the
general idea is the coupling of several C++-modules via Python, where each module does a different 
part of the alignment process. Good module combinations are finally coupled
under the roof of a single C++ application.
## Algorithmic Approaches
Many aligners conceptually incorporate the following 3 stage approach:
1. **Seeding**. Computation of perfect matches (or with few mismatches) using some form of indexing mechanism or hashing.
2. **Coupling**. Identification of promising seed subsets that could lead to an optimal alignment. A 
popular technique for this task is *chaining*.
3. **Dynamic Programming**. Closing of remaining gaps between seeds and extension of the outer 
end points given by the seeds.

MA stays with this basic pattern, but it relies on novel algorithmic solutions for the stages 1 and 2. 
MA introduces a divide and conquer approach for seeding on the foundation of the FMD-Index. The 
advantage of this variant of seeding is the reduction of the overall number of seeds compared to the
classical FMD-index based extension used in BWA-MEM. One could expect 
that such a reduction of the number of seeds hurts an aligner's accuracy. But this 
is only true, if we loose significant and insignificant seeds at equal rates. And, as this latter 
statement already suggests, the divide and conquer variant promotes significant seeds.\
In stage 2 we come up with two new ideas. First, the *strip of consideration* (SoC) and, 
second, *seed harmonization*. Both techniques are simple line-sweeps that iterate over a 
set of seeds. The SoC is used in order to quickly identify promising regions on the reference, 
while the seed harmonization purges contradicting 
seeds by relying on a virtual guide line. The highlights of these approaches are: No specially
tailored data structures required. Highly efficient. Easy to implement (roughly 50 lines Pseudocode).

## Python Integration

> **Python Integration is an optional feature of MA.** By default MA is built without Python support.

The Python integration of MA is done via [Boost.Python](https://www.boost.org/ "Boost.Python"). 
So, it is necessary to have a boost deployment with Python3 support. If there is Boost.Python
support, then `make` can be called with `WITH_PYTHON=1`. The folder *MA* of the repository 
contains a Python3 module `MA` that incorporates the aligner. The Python3 module relies on
the shared library *libMA.so*. Alignments using the Python integration 
can be done without a computational penalty, since Python is only responsible for an initial C++-module 
coupling, while all actual computations are done within the C++-modules. The idea is similar 
to the one used in the context of [TensorFlow](https://www.tensorflow.org "TensorFlow").

## Getting Started
MA is currently still in an evaluation and testing process. Feedback about bugs is highly welcomed. 
The below procedure was checked on a fresh *Debian 9.4 Stretch* with the following additional 
packages installed: `git`, `make` and `build-essential`.
### Installation
 
 Get the github clone and call make. (Available make switches are documented below):

    git clone https://github.com/ITBE-Lab/MA
    cd MA
    make

Possible Makefile switches:

- WITH_PYTHON=1\
    MA compiles as a shared library, which can be imported in Python.
    The MA executable is dynamically linked to the shared library, 
    so the folder comprising the shared library must be added to LD_LIBRARY_PATH.
    This options requires that $(BOOST_ROOT), $(PYTHON_INCLUDE) and $(PYTHON_LIB) are correctly set.
    *Tested using boost 1.65.1 (with the Boost.Python library) and Python3.*

- DEBUG=1\
    Compiles MA without optimizations and with assertions enabled. Further, there are multiple self-checks during runtime. The debug mode is primarily for bug hunting.

### First Steps ...
Test your installation with:

    ./ma -h

You should see the help screen of MA.


In order to perform alignments you need to **precompute an index** of your reference genome first.
This is done using following command:

    ./ma --genIndex -i <filename_of_genome> -x <filename_of_index>

`<filename_of_genome>` is the path to a FASTA file containing the reference genome and 
`<filename_of_index>` is the prefix that shall be used to store the index.

Now you can **align** using the following command:

    ./ma -x <filename_of_index> -i <fasta_in> -o <sam_out>

`<fasta_in>` is a FASTA or FASTQ file containing the queries.
`<sam_out>` is the filename of the output file that shall be created.
You can switch between **MA accurate** and **MA fast** by using the `-m acc` or 
`-m fast` switch, respectively.

## MA options

```

General options:
    -h, --help                     Display the complete help screen
        --genIndex                 Do FMD-index Generation. The -i and -x options specify the FASTA
                                   file used for index generation and the index prefix, respectively.
                                   If this option is not set, the aligner performs alignments. 

Necessary arguments for alignments:
    -x, --idx <prefix>             FMD-index used for alignments
    -i, --in <fname>               FASTA or FASTAQ input files.

Alignment options:
    -o, --out <fname>              Filename used for SAM file output. Default output stream is
                                   standard output.
    -t, --threads <num>            Use <num> threads. On startup MA checks the hardware and chooses 
                                   this value accordingly.
    -m, --mode [fast/acc]          Set operation modus for MA. 
                                   Default is 'fast'.
    -d, --noDP                     Switch that disables the final Dynamic Programming.
    -n, --reportN <num>            Report up to <num> alignments; 0 means unlimited.
                                   Default is 1.
    -s, --seedSet [SMEMs/maxSpan]  Selects between the two seeding strategies super maximal extend matches
                                   'SMEMs' and maximally spanning seeds 'maxSpan'. 
                                   Default is 'maxSpan'.
    -l, --minLen <num>             Seeds must have a minimum length of <num> nucleotides.
                                   Default is 16.
        --Match <num>              Sets the match score to <num>; <num> > 0.
                                   Default is 3. 
        --MisMatch <num>           Sets the mismatch penalty to <num>; <num> > 0.
                                   Default is 4.
        --Gap <num>                Sets the costs for opening a gap to <num>; <num> >= 0.
                                   Default is 6.
        --Extend <num>             Sets the costs for extending a gap to <num>; <num> > 0.
                                   Default is 1

Paired Read options:
    -p, --paUni                    Enable paired alignments and model the distance as uniform distribution.
                                   If set --in shall be used as follows: --in '<fname1>, <fname2>'.
    -P, --paNorm                   Enable paired alignment and Model the distance as normal distribution.
                                   If set --in shall be used as follows: '--in <fname1>, <fname2>'.
        --paIsolate <num>          Penalty for an unpaired read pair.
                                   Default is 17.
        --paMean <num>             Mean gap distance between read pairs.
                                   Default is 400.
        --paStd <num>              Standard deviation of gap distance between read pairs.
                                   Default is 150.

Advanced options:
        --giveUp <val>             Threshold with 0 <= <val> <= 1 used as give-up criteria.
                                   SoC's with accumulative seed length smaller than 
                                   'query_len * <val>' will be ignored.
                                   Reducing this parameter will decrease runtime but allows
                                   the aligner to discover more dissimilar matches.
                                   Increasing this parameter will increase runtime but might cause
                                   the aligner to miss the correct reference location.
                                   Default is 0.002.
        --maxTries <num>           Maximally <num> many SoC's are evaluated for a single alignment.
                                   Generally the best alignment is found in the best scored SoC.
                                   However, if the best alignment is from a very repetitive region,
                                   we might have to inspect several SoC's to find the optimal one.
                                   Default is 50.
        --minRefSize <num>         If the reference is smaller than <num> nt we disable all heuristics.
                                   Default is 10000000.
```

## Thanks ...

MA relies on the hard work of other projects. These are:

- Dynamic programming with adaptive band: https://github.com/ocxtal/libgaba
- Dynamic programming with static band: https://github.com/lh3/ksw2
- Cmd line option parser: https://github.com/jarro2783/cxxopts
- Parts of the code for FMD-index creation and extension were picked from https://github.com/lh3/bwa


Many thanks to the creators and contributors of the above projects ...

## Authors

MA was initiated and is maintained by Markus Schmidt and Arne Kutzner.
