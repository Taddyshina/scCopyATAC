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
  3. Go to scCopyATAC/scripts directory, dirictly start CNV analysis.```
  
  *** Type ```$ hint``` to test if HiNT successfully installed

### Download reference files used in scCopyATAC

1. Download hg19 reference [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html).  Download hg19 blacklists [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html).
2. Download hg38 reference [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html).  Download hg38 blacklists [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html).
3. Download mm10 reference [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html).  Download mm10 blacklists [HERE](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html).
4. Now only hg19, hg38 and mm10 available for scCopyATAC.

## Quick Start

* Download the test datasets from [HERE]()

### HiNT-PRE
HiNT pre: Preprocessing Hi-C data. HiNT pre does alignment, contact matrix creation and normalization in one command line.

```$ hint pre -d /path/to/hic_1.fastq.gz,/path/to/hic_2.fastq.gz -i /path/to/bwaIndex/hg19/hg19.fa --refdir /path/to/refData/hg19 --informat fastq --outformat cooler -g hg19 -n test -o /path/to/outputdir --pairtoolspath /path/to/pairtools --samtoolspath /path/to/samtools --coolerpath /path/to/cooler```

```$ hint pre -d /path/to/test.bam --refdir /path/to/refData/hg19 --informat bam --outformat juicer -g hg19 -n test -o /path/to/outputdir --pairtoolspath /path/to/pairtools --samtoolspath /path/to/samtools --juicerpath /path/to/juicer_tools.1.8.9_jcuda.0.8.jar```

use ```$ which samtools ``` ```$ which pairtools ``` ```$ which cooler ``` to get the absolute path of these tools, and ```/path/to/juicer_tools.1.8.9_jcuda.0.8.jar``` should be the path where you store this file

see details and more options

```$ hint pre -h ```

### HiNT-CNV
HiNT cnv: prediction of copy number information, as well as segmentation from Hi-C.

```$ hint cnv -m contactMatrix.cool -f cooler --refdir /path/to/refDir/hg19 -r 50 -g hg19 -n test -o /path/to/outputDir --bicseq /path/to/BICseq2-seg_v0.7.3 -e MboI```

```$ hint cnv -m /path/to/4DNFIS6HAUPP.mcool::/resolutions/50000 -f cooler --refdir /path/to/refDir/hg38 -r 50 -g hg38 -n HepG2 --bicseq /path/to/BICseq2-seg_v0.7.3 -e DpnII --maptrack 36mer```

```$ hint cnv -m /path/to/4DNFICSTCJQZ.hic -f juicer --refdir /path/to/refDir/hg38 -r 50 -g hg38 -n HepG2 --bicseq /path/to/BICseq2-seg_v0.7.3 -e DpnII```

```$ hint cnv -m /path/to/4DNFICSTCJQZ.hic -f juicer --refdir /path/to/refDir/hg38 -r 50 -g hg38 -n HepG2 --bicseq /path/to/BICseq2-seg_v0.7.3 -e DpnII --doiter```

```/path/to/BICseq2-seg_v0.7.3``` should be the path where you store this package

see details and more options

```$ hint cnv -h ```

### HiNT-TL
HiNT tl: interchromosomal translocations and breakpoints detection from
Hi-C inter-chromosomal interaction matrices.

```$ hint tl -m /path/to/data_1Mb.cool,/path/to/data_100kb.cool --chimeric /path/to/test_chimeric.sorted.pairsam.gz --refdir /path/to/refDir/hg19 --backdir /path/to/backgroundMatrices/hg19 --ppath /path/to/pairix -f cooler -g hg19 -n test -o /path/to/outputDir```

```$ hint tl -m /path/to/4DNFIS6HAUPP.mcool::/resolutions/1000000,/path/to/4DNFIS6HAUPP.mcool::/resolutions/100000 -f cooler --refdir /path/to/refDir/hg38 --backdir /path/to/backgroundMatrices/hg38 -g hg38 -n 4DNFICSTCJQZ -c 0.05 --ppath /path/to/pairix -p 12```

```$ hint tl -m /path/to/4DNFICSTCJQZ.hic -f juicer --refdir /path/to/refData/hg38 --backdir /path/to/backgroundMatrices/hg38 -g hg38 -n 4DNFICSTCJQZ -c 0.05 --ppath /path/to/pairix -p 12 -o HiNTtransl_juicerOUTPUT```

use ```$ which pairix ``` to get the absolute path of pairix

see details and more options

```$ hint tl -h ```

## Output of HiNT
### HiNT-PRE output
In the HiNT-PRE output directory, you will find

1. ```jobname.bam``` aligned lossless file in bam format
2. ```jobname_merged_valid.pairs.gz``` reads pairs in pair format
3. ```jobname_chimeric.sorted.pairsam.gz``` ambiguous chimeric read pairs used for breakpoint detection in [pairsam](https://github.com/mirnylab/pairtools) format
4. ```jobname_valid.sorted.deduped.pairsam.gz``` valid read pairs used for Hi-C contact matrix creation in [pairsam](https://github.com/mirnylab/pairtools) format
5. ```jobname.mcool``` Hi-C contact matrix in [cool](https://github.com/mirnylab/cooler) format
6. ```jobname.hic``` Hi-C contact matrix in [hic](https://github.com/aidenlab/juicer) format

### HiNT-CNV output
In the HiNT-CNV output directory, you will find

1. ```jobname_GAMPoisson.pdf``` the GAM regression result
2. ```segmentation/jobname_bicsq_allchroms.txt``` CNV segments with log2 copy ratio and p-values in txt file
3. ```segmentation/jobname_resolution_CNV_segments.png``` figure to visualize CNV segments
4. ```segmentation/jobname_bicseq_allchroms.l2r.pdf``` figure to visualize log2 copy ration in each bin (bin size = resolution you set)
5. ```segmentation/other_files``` intermediate files used to run BIC-seq
6. ```jonname_dataForRegression/*``` data used for regression as well as residuals after removing Hi-C biases

### HiNT-TL output
In the HiNT-TL output directory, you will find

1. ```jobname_Translocation_IntegratedBP.txt``` the final integrated translocation breakpoint
2. ```jobname_chrompairs_rankProduct.txt``` rank product predicted potential translocated chromosome pairs
3. ```otherFolders``` intermediate files used to identify the translocation breakpoints
