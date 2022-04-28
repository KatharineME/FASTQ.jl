# Benchmark

Fastq.jl has been benchmarked in the same fashion as entries into the [2016 PrecisionFDA Truth Challenge](https://precision.fda.gov/challenges/truth/results) and is therefore comparable to those entries.

The [2020 PrecisionFDA Truth Challenge V2](https://precision.fda.gov/challenges/10) was focused on difficult-to-map regions of GRCH38 including segmental duplications and MHC regions. In the future Fastq.jl will be benchmarked according to this challenge as well.

## PrecisionFDA Truth Challenge Overview

- Sequence analysis pipeline challenge where entires were competing on SNP and Indel recall, performance, and precision
- There were 35 entries
- Each challenge entry consisted of two VCFs: one for HG001 and one for HG002
- However only HG002 was used to evaluate the entries and choose the winners
- Genome in a Bottle (GiaB) consortium provided the truth data
- Global Alliance for Genomics and Health (GA4GH) provided software and best practices for comparisons
- Version 3.2.2 of HG002 truth data was used for evaluation
- The sex chromosomes did not have truth data, only chromosomes 1-22
- When evaluating entry VCFs they removed offending VCF lines, such as lines with "nan" in the REF column or non-diploid genotypes (0/1/2)

## Data

Data sources for the Fastq.jl benchmark

[HG002 raw data](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/) referenced from [here](https://github.com/genome-in-a-bottle/giab_data_indexes)

The folders under “HG002_HiSeq300x_fastq” each contain 20-30X sequencing (a single flow cell) and contain folders with fast files from each library, which can be combined for most purposes. Samples A-L are the six vials of starting material. Each sample has two technical replicates, hence Sample A1 and Sample A2.

[HG002 Truth VCF and BED](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/)

## Software

GA4GH was responsible for defining best practices and software specifications for the PrecisionFDA Truth Challenge. GA4GH created a prototype benchmarking workflow that consisted of specific implementations of RTG's vcfeval tool and Illumina's hap.py tool.

RTG’s vcfeval was used for VCF comparison. vcfeval generates an intermediate VCF which is further quantified by hap.py (the HAP-207 version, with the engine set to a GA4GH-specific version of vcfeval). hap.py’s quantity tool counts and stratifies variants by type.

This benchmark workflow is documented [here](https://github.com/ga4gh/benchmarking-tools/tree/master/doc/ref-impl).

![pipeline](media/benchmark_pipeline.png)
￼

- The intermediate.vcf created by vcfeval must have two columns names: TRUTH and QUERY with these FORMAT annotations: `##FORMAT=<ID=BK,Number=1,Type=String,Description="Sub-type for decision (match/mismatch type)">`
- The stratification bed files tell hap.py which variant is what type. hapy.py uses that information to give performance metrics in each variant type category
  - The new stratification bed files according to ga4gh: https://github.com/genome-in-a-bottle/genome-stratifications
- How you would call hap.py to evaluate a query.vcf

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz.tbi
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_v4.1_SmallVariantDraftBenchmark_12182019/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.bed

hap.py HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz query.vcf.gz -f HG002_GRCh37_1_22_v4.1_draft_benchmark.bed -o benchmarking-output

Vcfeval

- Before running need to convert reference faster to sdf format like this:
  `rtg format -o hg19.sdf hg19.fa` - that command creates a folder hg38.sdf which must be passed to hap.py `--engine-vcfeval-template hg38.sdf`

Hap.py -https://github.com/Illumina/hap.py

This is command to use with the exception that rtg wasn’t installed via hap.py:

`hap.py truth.vcf.gz query.vcf.gz -f conf.bed.gz -o ./test -V --engine-vcfeval-path /path/to/rtg --engine-vcfeval-template /path/to/hg38.sdf`

- `-f` is for passing a bed file of confident call regions, hap.py will be able to say that variant calls different from these are false positives
- `-o` specifies an output file prefix
- `-V` writes an annotated VCF
- `--engine-vcfeval-path` specifies path to rtg installation
- `--engine-vcfeval-template` specifies path to SDF used for vcfeval run

## Metrics

Highest SNP performance - highest SNP F-score
Highest SNP recall
Highest SNP precision
Highest Indel performance - highest Indel F-score
Highest Indel recall
Highest Indel precision

![metrics](media/precisionfda_metrics.png)

- True Positive / TP: present in both truth and query
- False Positive / FP: present only in the query
- False Negative / FN: present only in the truth
- Not-assessed / N: call was not assigned a match status
