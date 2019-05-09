# <img src="https://raw.githubusercontent.com/ITBE-Lab/MA/release/MA.png" align="center" width="90"> The Modular Aligner
MA is a novel open source application for the efficient and accurate alignment of short and long reads of various sequencers. The approach has a highly modular architecture and everyone is invited to propose/integrate new modules. 
The design aims at a smooth Python integration, while keeping the performance delivered by C++. So, the
general idea is the coupling of several C++-modules via Python, where each module does a different 
part of the alignment process. Good module combinations are finally coupled
under the roof of a single C++ application.

## Algorithmic Approach
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


## Getting Started

### Installer for MA on Windows 7/10 (64 bit)
An **installer for MA on Windows 7/10 (64 bit)** is available
[**here**](http://itbe.hanyang.ac.kr/software/ma/ "MA for Windows 7/10 (64 bit)").\
After installing MA, start the GUI and build an FMD-index using the menu entry 
`Genome->Create Index` or the `F2` key.
Then, follow the `Quickstart` description in the upper-right corner for computing alignments.

### Installation on Linux (Debian and Ubuntu)
The below procedure was checked on a fresh *Debian 9.4 Stretch* with the following additional 
packages installed: `git`, `cmake`, `make` and `build-essential`.\
Get the github clone and use cmake (at least version 3.8) for building. This can be done via: 

    git clone https://github.com/ITBE-Lab/MA
    cd MA
    mkdir build
    cd build
    cmake ..
    make

### Bioconda

MA can be installed using bioconda via:

    conda install ma

### Compilation on Windows 7/10
Get the github clone via:

    git clone https://github.com/ITBE-Lab/MA

Then, MA can be build via the `cmake` support of Microsoft Visual Studio 2017.


## First Steps with the command line version ...
Test your installation with:

    ./maCMD -h

You should see the help screen of MA.

In order to perform alignments you need to **precompute an index** of your reference genome first.
This is done using following command:

    ./maCMD --Create_Index <fasta_file_name>,<output_folder>,<index_name>

`<fasta_file_name>`: filepath (filename) of the Fasta file holding the genome used for index creation. \
`<output_folder>`: folderpath (foldername) of the location used for index storage. \
`<index_name>`: name used for identifying the new FMD-Index. In the context of alignments,
the genome-name is used for FMD-index selection.

Now you can **align** using the following command:

    ./maCMD -x <index_name> -i <fasta_in> -o <sam_out>

`<fasta_in>` is a FASTA or FASTQ file containing the queries.
`<sam_out>` is the filename of the output file that shall be created.


## Optional Features

### Python Integration
MA is built with Python support, if cmake can locate a valid python installation (at least version 3.5).
The Python integration of MA is done via [pybind11](https://github.com/pybind/pybind11 "pybind11"). 
Alignments using the Python integration can be done without a computational penalty, 
since Python is only responsible for an initial C++-module 
coupling, while all actual computations are done within the C++-modules. The idea is similar 
to the one used in the context of [TensorFlow](https://www.tensorflow.org "TensorFlow").

### zLib Support
For getting zLib support cmake has to find the appropriate zLib libraries.
On Debian/Ubuntu you have to run `apt-get install zlib1g-dev`.
On Windows the libraries are available via the [zLib homepage](https://www.zlib.net/ "zLib homepage").

### Graphical User Interface (GUI) Support via wxWidgets
The GUI version is additionally build to the comand line version,
if `cmake` can find support for [wxWidgets](https://www.wxwidgets.org/ "wxWidgets Homepage").
On Debian/Ubuntu, wxWidgets can be installed via `apt-get install libwxgtk3.0-gtk3-dev`.
On Windows wxWidgets must be downloaded from the [wxWidgets homepage](https://www.wxwidgets.org/ "wxWidgets Homepage").


## Changelog:

Version 1.1.1:

* MA can now recognize short inversions too. A z-drop based heuristic is used as algorithmic 
approach for inversion recognition. 
Regions between seeds showing a z-drop are reverse complemented and realigned using DP. 
If the DP-score is sufficiently high, the region is accepted as inversion and a supplementary alignment 
is provided using the computed DP-path. Similar approaches are used by other aligners too.
* Minor bug fixes.


## MA Options

```
Available presettings:
    -p, --Presetting <name> [Default]              Optimize aligner parameters for a selected
                                                   sequencing technique. Available presettings are:
                                                   'Default', 'Illumina', 'Illumina_Paired',
                                                   'Nanopore', 'PacBio'.

General options: (these options are not affected by presettings)
    -x, --Index <file_name> []                     Filename of FMD-index. (A FMD-index can be
                                                   generated via the --Create_Index option.) This
                                                   option must be set.

    -i, --In <file_name> []                        Filenames of Fasta/Fastq files containing reads.
                                                   gz-compressed files are automatically decompressed.
                                                   Multiple files can be specified by a comma
                                                   separated list. At least one file name must be
                                                   provided.

    -m, --Mate_In <file_name> []                   Filenames of the mates in the case of paired reads.
                                                   If this option is set, the aligner switches to
                                                   paired mode automatically. The number of reads
                                                   given as mates must match the accumulated number of
                                                   reads provided via the 'in'-option.

    -X, --Create_Index <fasta_file_name,output_folder,index_name> []
                                                   Generate a FMD-index for a Fasta file.
                                                   'fasta_file_name' has to be the file-path of the
                                                   Fasta file holding the genome used for index
                                                   creation. 'output_folder' is the folder-path of the
                                                   location used for index storage. 'index_name' is
                                                   the name used for identifying the new FMD-Index. In
                                                   the context of alignments, the genome-name is used
                                                   for FMD-index selection.

    -o, --SAM_File_name <file_name> [ma_out.sam]
                                                   Name of the SAM file alignments shall be written
                                                   to.

    -t, --Number_of_Threads <int> [32]             Number of threads used in the context of
                                                   alignments. This options is only available, if 'use
                                                   all processor cores' is off.

    -h, --Help <bool> [true]                       Print the complete help text.

Paired Reads options:
    --Use_Paired_Reads <bool> [false]              If your reads occur as paired reads, activate this
                                                   flag.

    -d, --Mean_Distance_of_Paired_Reads <double> [400]
                                                   Two reads can be paired if they are within mean +-
                                                   (standard deviation)*3 distance from one another on
                                                   the expected strands (depends on Use Mate Pair
                                                   on/off) Used in the context of the computation of
                                                   the mapping quality and for picking optimal
                                                   alignment pairs.

    -S, --Standard_Deviation_of_Paired_Reads <double> [150]
                                                   <val> represents the standard deviation for the
                                                   distance between paired reads. Used in the context
                                                   of the computation of the mapping quality and for
                                                   picking optimal alignment pairs.

    --Score_Factor_for_Paired_Reads <double> [1.25]
                                                   This factor is multiplied to the score of
                                                   successfully paired reads. Used in the context of
                                                   the computation of the mapping quality and for
                                                   picking optimal alignment pairs. [val] < 1 results
                                                   in penalty; [val] > 1 results in bonus.

    --Check_for_Consistency <bool> [false]         Check if both paired read files comprise the same
                                                   number of reads. (Intended for debugging.)

Seeding options:
    -s, --Seeding_Technique <name> [maxSpan]       Technique used for the initial seeding. Available
                                                   techniques are: maxSpan and SMEMs.

    -l, --Minimal_Seed_Length <int> [16]           All seeds with size smaller than 'minimal seed
                                                   length' are discarded.

    --Minimal_Ambiguity <int> [0]                  During the extension of seeds using the FMD-index:
                                                   With increasing extension width, the number of
                                                   occurrences of corresponding seeds on the reference
                                                   montonically decreases. Keep extending, while the
                                                   number of occurrences is higher than 'Minimal
                                                   Ambiguity'. (For details see the MA-Handbook.)

    --Maximal_Ambiguity <int> [100]                Discard seeds that occur more than 'Maximal
                                                   ambiguity' time on the reference. Set to zero to
                                                   disable.

    --Skip_Ambiguous_Seeds <bool> [false]          Enabled: Discard all seeds that are more ambiguous
                                                   than [Maximal Ambiguity]. Disabled: sample [Maximal
                                                   Ambiguity] random seeds from too ambiguous seeds.

    --Seeding_Drop-off_A_-_Minimal_Seed_Size <int> [15]
                                                   Heuristic runtime optimization: For a given read R,
                                                   let N be the number of seeds of size >= [val].
                                                   Discard R, if N < [length(R)] * [Seeding drop-off
                                                   B].

    --Seeding_Drop-off_B_-_Factor <double> [0.005]
                                                   Heuristic runtime optimization: Factor for seed
                                                   drop-off calculation. For more information see
                                                   parameter [Seeding drop-off A].

Strip of Consideration options:
    -N, --Maximal_Number_of_SoC's <int> [30]       Only consider the <val> best scored SoC's. 0 = no
                                                   limit.

    -M, --Minimal_Number_of_SoC's <int> [1]        Always consider the first <val> SoC's no matter the
                                                   Heuristic optimizations. Upping this parameter
                                                   might improve the output of supplementary
                                                   alignments and therefore successive SV calling.

    --Fixed_SoC_Width <int> [0]                    Set the SoC width to a fixed value. 0 = use the
                                                   formula given in the paper. This parameter is
                                                   intended for debugging purposes.

SAM Output options:
    -n, --Maximal_Number_of_Reported_Alignments <int> [0]
                                                   Do not output more than <val> alignments. Set to
                                                   zero for unlimited output.

    --Minimal_Alignment_Score <int> [75]           Suppress the output of alignments with a score
                                                   below val.

    --Omit_Secondary_Alignments <bool> [false]     Suppress the output of secondary alignments.

    --Omit_Supplementary_Alignments <bool> [false]
                                                   Suppress the output of supplementary alignments.

    --Maximal_Supplementary_Overlap <double> [0.1]
                                                   An non-primary alignment A is considered
                                                   supplementary, if less than val percent of A
                                                   overlap with the primary alignment on the query.
                                                   Otherwise A is considered secondary.

    --Number_Supplementary_Alignments <int> [1]
                                                   Maximal Number of supplementary alignments per
                                                   primary alignment.

    --Emulate_NGMLR's_tag_output <bool> [false]
                                                   Output SAM tags as NGMLR would. Activate this if
                                                   you want to use MA in combination with Sniffles.
                                                   Enableing this will drastically increase the size
                                                   of the output file.

    --Output_long_cigars_in_CG_tag <bool> [true]
                                                   Some software crashes if a cigar is too long.
                                                   Enabeling this flag makes MA output that the entire
                                                   read was soft clipped in the regular cigar field if
                                                   the cigar would exceed 65536 operations. The actual
                                                   cigar is then given in the CG:B:I tag as a comma
                                                   seperated binary list.

Heuristics options:
    --SoC_Score_Drop-off <double> [0.1]            Let x be the maximal encountered SoC score. Stop
                                                   harmonizing SoC's if there is a SoC with a score
                                                   lower than <val>*x.

    --Minimal_Harmonization_Score <int> [18]       Discard all harmonized SoC's with scores lower than
                                                   <val>. Only keep detected inversions with a score >=
                                                   <val> * [Match Score].

    --Relative_Minimal_Harmonization_Score <double> [0.002]
                                                   Discard all harmonized SoC's with scores lower than
                                                   length(read)*<val>.

    --Harmonization_Drop-off_A_-_Score_Difference <double> [0.0001]
                                                   Let x be the maximal encountered harmonization
                                                   score. Stop harmonizing further SoC's if there are
                                                   <Harmonization Drop-off B> SoC's with lower scores
                                                   than x-<readlength>*<val> in a row.

    --Harmonization_Drop-off_B_-_Lookahead <int> [3]
                                                   See Harmonization Drop-off A.

    --Harmonization_Score_Drop-off_-_Minimal_Query_Length <int> [800]
                                                   For reads of length >= [val]: Ignore all SoC's with
                                                   harmonization scores lower than the current maximal
                                                   score. 0 = disabled.

    --Artifact_Filter_A_-_Maximal_Delta_Distance <double> [0.1]
                                                   Filter seeds if the difference between the delta
                                                   distance to it's predecessor and successor is less
                                                   then [val] percent (set to 1 to disable filter) and
                                                   the delta distance to it's pre- and successor is
                                                   more than [Artifact Filter B] nt.

    --Artifact_Filter_B_-_Minimal_Delta_Distance <int> [16]
                                                   See Artifact Filter A

    --Pick_Local_Seed_Set_A_-_Enabled <bool> [false]
                                                   <val> = true enables local seed set computation.

    --Pick_Local_Seed_Set_B_-_Optimistic_Gap_Estimation <bool> [true]
                                                   After the harmonization MA checks weather it is
                                                   possible to compute a positively scored alignment
                                                   from the seed set. Gaps between seeds can be
                                                   estimated in two ways: Optimistic [true]: Assume
                                                   that the gap can be filled using merely matches and
                                                   a single insertion/deletion. Pessimistic [false]:
                                                   Assume that the gap can be filled using matches and
                                                   mismatches that add up to a score of 0 and a single
                                                   insertion/deletion.

    --Pick_Local_Seed_Set_C_-_Maximal_Gap_Penalty <int> [100]
                                                   Maximal Gap cost penalty during local seed set
                                                   computation.

    --Maximal_Gap_Size <int> [20]                  If the gap between seeds is larger than <val> on
                                                   query or reference, the dual extension process is
                                                   used to fill the gap. Dual extension is more
                                                   expensive if the extension does not Z-drop, but
                                                   more efficient otherwise.

    --Minimum_Genome_Size_for_Heuristics <int> [10000000]
                                                   Some heuristics can only be applied on long enough
                                                   genomes. Disables: SoC score Drop-off if the genome
                                                   is shorter than <val>.

    --Disable_All_Heuristics <bool> [false]        Disables all runtime heuristics. (Intended for
                                                   debugging.)

Dynamic Programming options:
    --Match_Score <int> [2]                        Match score. (Used in the context of Dynamic
                                                   Programming and for SoC width computation.)

    --Mismatch_Penalty <int> [4]                   Penalty for mismatch.

    --Gap_penalty <int> [4]                        First penalty for gap opening. (Two piece affine
                                                   gap costs)

    --Extend_Penalty <int> [2]                     First penalty for gap extension. (Two piece affine
                                                   gap costs)

    --Second_Gap_Penalty <int> [24]                Second penalty for gap opening. (Two piece affine
                                                   gap costs)

    --Second_Extend_Penalty <int> [1]              Second penalty for gap extension. (Two piece affine
                                                   gap costs)

    --Padding <int> [1000]                         If an alignment does not reach its read's
                                                   endpoints, the missing parts can be computed via
                                                   dynamic programming. If the length of the missing
                                                   parts is smaller than 'Padding', dynamic
                                                   programming is used to extend the alignment towards
                                                   the endpoints of the read. Otherwise, the unaligned
                                                   parts of the read are ignored and the alignment
                                                   stays unextended.

    --Bandwidth_for_Extensions <int> [512]         Bandwidth used in the context of extending an
                                                   alignment towards the endpoints of its read. (See
                                                   'Padding')

    --Minimal_Bandwidth_in_Gaps <int> [20]         Gaps between seeds are generally filled using
                                                   dynamic programming. This option determines the
                                                   minimal bandwidth used in the context of bridging
                                                   gaps. More details can be found in the MA-Handbook.

    --Z_Drop <int> [200]                           If the running score during dynamic programming
                                                   drops faster than <val> stop the extension process.

    --Detect_Small_Inversions <bool> [false]       Use DP to search for small inversions that do not
                                                   contain any seeds.

    --Z_Drop_Inversions <int> [100]                Check for an inversion if the running score during
                                                   dynamic programming drops faster than <val>.

```

## Thanks ...

MA relies on the hard work of other projects. These are:

- Dynamic programming with static band: https://github.com/lh3/ksw2
- Parts of the code for FMD-index creation and extension were picked from https://github.com/lh3/bwa

Many thanks to the creators and contributors of the above projects ...

## Citing MA

For citing MA, please use:

Schmidt, M., Heese, K. & Kutzner, A. Accurate high throughput alignment via line sweep-based seed processing. Nature Communications 10, 1939, doi:10.1038/s41467-019-09977-2 (2019).

## Authors

MA was initiated and is maintained by Markus Schmidt and Arne Kutzner.
