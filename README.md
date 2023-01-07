# scCopyATAC
## A computational method for delineating copy number and clonal substructure from tumor single-cell ATAC-seq data without normal reference

## Summary
**scCopyATAC**, a computational method to detect CNVs and delineat clonal substructure from single cell ATAC sequencing data without normal reference. scCopyATAC has two main steps: **scCopyATAC-PRE**, **scCopyATAC-CNV**. scCopyATAC-PRE preprocesses single cell ATAC sequencing data and computes the raw CNV signals matrix. And do clustering for single cell by chromotin accessibility, which is a paperation for next step detecting diploid baseline. Then scCopyATAC-CNV detect potential diploid baseline by gaussian mixture model and redefine raw CNV signals get the final CNV value. scCopyATAC using stand cellranger output as input.

#### Overview of scCopyATAC workflow: 
<img src="https://github.com/Taddyshina/scCopyATAC/blob/main/figures/scCopyATAC_pipeline.png" width="800">

## Installation

### Dependencies
R and R packages

1. [R >= 3.4](https://www.r-project.org/)
2. [ChIPpeakAnno](https://www.bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html), [Signac](https://cloud.r-project.org/web/packages/Signac/index.html), [Seurat](https://cloud.r-project.org/web/packages/Seurat/index.html), [GenomeInfoDb](https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html), [patchwork](https://cran.r-project.org/web/packages/patchwork/index.html), [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), [matrixStats](https://cran.r-project.org/web/packages/matrixStats/index.html), [readr](https://cran.r-project.org/web/packages/readr/index.html), [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), [magrittr](https://cran.r-project.org/web/packages/magrittr/index.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).

### Install scCopyATAC

* scCopyATAC is a install free pipeline which can be dirictly used. 
  1. Install scCopyATAC dependencies
  2. Download scCopyATAC ```git clone https://github.com/Taddyshina/scCopyATAC.git```
  3. Go to scCopyATAC/scripts directory, dirictly start CNV analysis.

### Download reference files used in scCopyATAC

1. Download hg19 reference [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html).  Download hg19 blacklists [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html).
2. Download hg38 reference [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html).  Download hg38 blacklists [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html).
3. Download mm10 reference [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html).  Download mm10 blacklists [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html).
4. Now only hg19, hg38 and mm10 available for scCopyATAC.

## Quick Start

* Download the test datasets from [HERE]()

### scCopyATAC-PRE
scCopyATAC-PRE: Preprocesses single cell ATAC sequencing data and computes the raw CNV signals matrix. And do clustering for single cell by chromotin accessibility in one command line.

```$ Rscript scCopyATAC-PRE.r 10X_h5 fragment metadata outputPrefix subBC species blacklists CGneighbors Win_size step_length```

**10X_h5**              :h5 file of filtered_peak_bc_matrix, cellRanger output.

**metadata_pathway**           :csv file of singlecell metadata, cellRanger output.

**fragment**			:fragments.tsv, cellRanger output.

**outPrefix**			:Prefix of output.It should contain directory path. (eg. /out/dir/namePrefix)

**subBC**				:Barcode list file. It should be one column and no header.

**species**             :Species and references version. "hg19" for hg19; "hg38" for hg38; "mm10" for mm10.

**blacklists**          :matched blacklist files of species.

**CGneighbors**         :Considered CG contents to caculate CNV.(default value is 100)

**Win_size**            :windows size to caculate CNV.(default value is 10e6)

**step_length**         :step length of windows(default value is 2e6)



### scCopyATAC-CNV
scCopyATAC-CNV: detect potential diploid baseline by gaussian mixture model and redefine raw CNV signals get the final CNV value

```$ Rscript scCopyATAC-CNV.r input sample resolution```

**input**              :RDS output of scCopyATAC-pre

**sample**           :sample type cell line or tumor.(eg. Cell_line, Tumor)

**resolution**			:Clustering resolution of subclone clustering.


## Output of scCopyATAC
### scCopyATAC-PRE output
In the scCopyATAC-PRE output directory, you will find

1. ```sample_metadata_S1.txt``` cell info contains cell barcodes, cell clusters base on chromotin accessibility.
2. ```sample_windows.txt``` windows info contains cnv windows chromosomo, start pos, end pos and CG content.
3. ```sample_S1.rds``` files for next step scCopyATAC-CNV.
4. ```sample_raw_CNV_heatmap.png``` heatmap of raw CNV signals.
5. ```sample_raw_CNV.mtx``` raw CNV cellxwindows mtx.

### scCopyATAC-CNV output
In the scCopyATAC-CNV output directory, you will find

1. ```sample_metadata_S2.txt``` cell info contains cell barcodes, cell clusters base on chromotin accessibility, predicated cell types and predicated subclone.
2. ```sample_final_CNV_heatmap.png``` heatmap of final CNV signals.
3. ```sample_final_CNV_line_profile.png``` line profile of final CNV signals.
4. ```sample_final_subclone_tree.png``` clonal relationship of subclone
5. ```sample_final_CNV.mtx``` final CNV cellxwindows mtx.
