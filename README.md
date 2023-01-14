## NGS tools
 

#### Command line `Rscript` tools:
  - **frag_coverage_script.R**: generates insertion size coverage bigwig from a paired-end bam file.
  - **ATAC_QC.R**: process a bam file and called peaks each time. Append reads quality (mapped rate) and enrichment (FRiP and TSS fraction) to `ATAC_report.txt`.


  
#### `R` functions and modules:
  - **HMM.R**: genome segmentation with HMM gaussian model. Written in `R` and `Rcpp` as an internal dependency for [`TU filter`](https://github.com/shaorray/TU_filter) alternative to `GenoSTAN`'s HMM utility.
