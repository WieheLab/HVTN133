# HVTN133

This repository contains analysis scripts used in the manuscript, Williams et al. "Vaccine Induction in Humans of Polyclonal HIV-1 Heterologous Neutralizing Antibodies".

## About

This repository makes available scripts that were used for two downstream analyses in the Williams et al. manuscript:
* Analysis of 10X Single Cell VDJ sequences
* Calculation of lipid insertion propensity scores for amino acid sequences from CDR regions

## Table of contents

> * [HVTN133](#hvtn133)
>   * [About](#about)
>   * [Table of contents](#table-of-contents)
>   * [Analysis of 10X Single Cell VDJ sequences](#analysis-of-10x-single-cell-vdj-sequences)
>   * [Lipid insertion propensity scores](#lipid-insertion-propensity-scores)
>   * [References](#references)

## Analysis of 10X Single Cell VDJ sequences


### Workflow:

1. CellRanger (v6.0.0) mkfastq was run on each library.
2. CellRanger vdj and count were run for the immune profiling and gene expression libraries, respectively. 
3. Immune Profiling (VDJ):

   a. All contigs from CellRanger were filtered for cell calls and high confidence. We also removed contigs with low UMI counts compared to the other contigs in that cell and for that chain type (all_contig_annots_umi_filter.R).

   b. Data was then run through Cloanalyst and functional and productive contigs with one heavy and one light chain were filtered.

   c. Clonal partitioning was then performed for each chain.

   d. The merge_clone_info.R script was then run to combine heavy and light chain immunogenetics and clonal assignments.

   e. Immunogenetic information of interest was then analyzed using the HVTN133_MPER_VDJ_Plots.R and VH7-4-1_Usage_Plot.R scripts.
4. Gene Expression: data was filtered, integrated, visualized, and analyzed using Seurat.

## Lipid insertion propensity scores

### Workflow

Hydropathy analysis online tool MPEx [[1]](#1) (Membrane Protein Explorer) v3.3.0 was used with the water-interface scale to calculate lipid insertion propensity scores as the sum of ΔGwif, the free energy of transfer of an amino acid from water to POPC interface [[2]](#2), over all amino acids in the CDR regions of both light and heavy chains. CDR amino acid positions were defined using the software ANARCI [[3]](#3) with the IMGT scheme. Custom perl and R scripts were used to compute scores for all CDRs individually.

### Usage

1. MPEx and ANARCI are run on sequences of interest. Note: MPEx is run in batch mode on multiple sequences.
2. Perl script parse_MPEx.pl is run on MPEx results to split the results into one file per sequence with the score for each individual amino acid.
3. R markdown script breakdown_dwif_by_CDR.Rmd is then run to combine the ANARCI results to split the sequences into CDR regions and calculate a combined lipid insertion propensity score for each region.

This workflow is adapted from a workflow described in a previous manuscript [[4]](#4).

## References
<a id="1">[1]</a> 
C. Snider, S. Jayasinghe, K. Hristova, S. H. White, MPEx: a tool for exploring membrane proteins. Protein Sci 18, 2624-2628 (2009).

<a id="2">[2]</a> 
W. C. Wimley, S. H. White, Experimentally determined hydrophobicity scale for proteins at membrane interfaces. Nat Struct Biol 3, 842-848 (1996).

<a id="3">[3]</a> 
J. Dunbar, C. M. Deane, ANARCI: antigen receptor numbering and receptor classification. Bioinformatics 32, 298-300 (2016).

<a id="4">[4]</a> 
R. Zhang, L. Verkoczy, K. Wiehe, S. Munir Alam, N. I. Nicely, S. Santra, T. Bradley, C. W. Pemble IV, J. Zhang, F. Gao, D. C. Montefiori, H. Bouton-Verville, G. Kelsoe, K. Larimore, P. D. Greenberg, R. Parks, A. Foulger, J. N. Peel, K. Luo, X. Lu, A. M. Trama, N. Vandergrift, G. D. Tomaras, T. B. Kepler, M. A. Moody, H. X. Liao, B. F. Haynes, Initiation of immune tolerance-controlled HIV gp41 neutralizing B cell lineages. Sci. Transl. Med. 8, 336ra362 (2016).