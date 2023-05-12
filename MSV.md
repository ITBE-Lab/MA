# MSV - The Modular SV Caller
MSV is the prototype of a SV caller that shall demonstrate the viability of the approaches proposed in the manuscript [[Preprint] State-of-the-art structural variant calling: What went conceptually wrong and how to fix it?](https://biorxiv.org/cgi/content/short/2021.01.12.426317v1 "bioRxiv preprint") Currently, this prototype is usable with small genomes but not ready for the use with large genomes as e.g. the human genome. Further, it does not deliver calls in the VCF format but stores computed adjacency matrix entries in a database instead. 

## Algorithmic Approach
[There is a blog about the algorithmic approach of MSV for nested SV.](http://itbe.hanyang.ac.kr/ak/MSV/ "MSV - blog")

MSV breaks with the state-of-the-art approaches for split-read and read-pair SV calling in two ways:
- Genomic rearrangements (SV) are detected via maximal extended matches (MEMs) instead of alignments. This approach eliminates the concealing effects that alignments have on SV.
- Genomic rearrangements are represented using adjacency matrices of a skew symmetric graph model for avoiding ambiguities that are inherent to the representation of SV via deletions, insertions and other forms of atomic SV. Further, our skew symmetric graph model unifies forward strand and reverse strand.

By breaking with existing approaches for SV calling, we are capable of unambiguously representing arbitrarily complex, nested SV. This is demonstrated in the manuscript by the lossless representation of one yeast genome (UFRJ50816) via a set of SV calls (i.e. an adjacency matrix) for another yeast genome (YPS138).

## Getting started
For installing MSV follow the commands given here.

## Running MSV
Currently, MSV is executed and controlled via python3. As backend, it relies on a comprehensive C++ library and a PostgreSQL database. Visualizations are done via the bokeh library for python. The below Python script shows a minimal example for getting MSV running on a new dataset. It executes three stages:
- Import reads into a database.
- Compute areas of adjacency matrix entries for all reads.
- Cluster overlapping areas.
<!-- Comment inserted to render code *outside* the list -->
    from MSV import *

    # select DB-Name
    db_name = "your-db-name"

    # load reference genome
    pack = Pack()
    pack.load("path-to-reference-genome-pack")
    # OR 
    # pack.append_fasta_file("path-to-reference-genome-fasta")
    # pack.store("path-to-reference-genome-pack)

    # select the parameter set
    parameter_set = "SV-PacBio" # OR parameter_set = "SV-Illumina"
    param = ParameterSetManager()
    param.set_selected(parameter_set)

    # load reads into DB
    seq_ids = [
                insert_reads_path_string_vec(param, 
                                            db_name,
                                            "read-dataset-name",
                                            ["read-fastq-or-fasta-file-paths"],
                                            ["mate-fastq-file-paths-or-none-for-long-reads"],
                                            coverage=100 # expected coverage
                                            )
                    ]

    # create Minimizer index
    mm_index = MinimizerIndex(param, pack.contigSeqs(), pack.contigNames())

    # compute entry areas of reads
    jump_id = compute_sv_jumps(param, mm_index, pack, db_name, seq_ids)

    # cluster entry areas for reads
    sv_caller_run_id = sweep_sv_jumps(param, db_name, jump_id, "run-name", "run-description", [0], pack)

After the execution of the above script, the clustered adjacency matrix entries are in the "sv_call_table" table of the "your-db-name" database. You can inspect the entries by using the PostgreSQL tool pgAdmin or by using our bokeh based visualization tool as described in the next section.

## Visualizing adjacency matrices
MSV represents SV calls (genomic rearrangements) via a skew-symmetric graph model. Adjacency matrices of this graph model can be visualized using a Python script that relies on bokeh as backend.

You can start the visualization tool from MA's main directory via the command:

    cd build/MSV/sv_visualization
    bokeh serve --show /bokeh_server.py
This opens a browser window.

Using the top-right drop-down buttons, select your dataset, run-id and ground-truth.
- Crosses represent adjacency matrix entries of the selected run.
- Circles represent ground-truth entries.
- You can zoom into the adjacency matrix by scrolling in the top-left plot; eventually the entry areas of reads and their seeds (leftmost and bottom-left plot) will render.
- You can click on seeds (horizontal or vertical lines in the bottom-left or top-left plot, respectively) for visualizing all seeds of a read in a diagrammatic dot-plot (center-right plot).
- By toggling the "Render true-positives / false-negatives / false-positives" buttons and using the “score” slider, you can filter down the dataset.
- Increasing the “max. render” slider will render more detail. (Here, large values will slow down visualization a lot.)

For computing recall & accuracy, adjust the Blur setting slider in the bottom right to
- 0: for Illumina reads
- 100: for PacBio reads

After adjusting the slider, click the button "Compute Stats", which is the fourth button above the slider. After some time, the curves appear in the plot window below the slider. Using the tabs on top of the plot window, you can switch between the raw numbers (Tab “Min Score”) and the recall & accuracy rates (Tab “Recall & Accuracy”). You can zoom into the plots using the mouse wheel. Note: You must be zoomed in enough to see individual adjacency matrix entry crosses and ground-truth circles for computing recall & accuracy.

## Recreating the experiments of the manuscript
The experiments of [[Preprint] State-of-the-art structural variant calling: What went conceptually wrong and how to fix it?](https://biorxiv.org/cgi/content/short/2021.01.12.426317v1 "bioRxiv preprint")
can be recreated using the [MSV-EVAL repository](https://github.com/ITBE-Lab/MSV-EVAL "MSV-EVAL").

## Changelog

Version 2.0.2:

* 2nd revision for MSV
* The confidence scores for variant calls now make use of coverage information. Before, the conficance of a call was the number of reads supporting the call. Now it is the number of reads supporting over the the coverage in the region of the call. This helps with making calls in real world data, where coverage is uneven.

Version 2.0.1:

* 1st Revision for MSV - only minor changes.