## Trim

### fastp

Trims and preprocesses read files.

- `detect_adapter_for_pe` Detect the adapter sequence. This is done by default with single-end sequencing, but for paired end, you must pass this option.
- `json` Specify path for json report
- `html` Specify path for html report

## Align

### minimap2

- `-x` A recommended meta flag that specifies alignment chain bandwidth, elongation, discard, how many secondary alignment should be output, and more.
- `-a` Generates CIGAR and output alignments in SAM format.
- `--sr` Enable short read mode.
- `--splice` Enable splice alignment mode.
- `-uf` Use transcript strand to find canonical splice sites.
- `-t NUM` Number of threads.
- `-K NUM` Number of bases loaded into memory to process in a mini-batch [500M]. A large NUM passed here helps with load balancing in the multi-threading mode.
- `-R` SAM read group line in a format like @RG\\tID:foo\\tSM:bar [].

### samtools fixmate

Corrects flaws in read-pairing that the aligner may have introduced. It ensures the SAM FLAG, RNEXT, PNEXT, and TLEN fields are correct and consistent.

- `-u` Uncompressed output.
- `-m` Add mate score tag.
- `-` A synonym for stdin and stdout in samtools.

### samtools sort

Orders aligned reads by chromosome and coordinate. `sort` is highly parallel so adding threads is beneficial.

- `-u` Uncompressed output.
- `-l 1` Fastest level of BAM compressed output.
- `-T PATH` Specify a temporary storage directory like `/tmp/example_prefix/`.

### samtools markdup

Marks duplicate alignments in a coordinate sorted file using flags added by `fixmate -m`. Duplicates are defined as aligned reads whose 5 prime coordinates and orientation (forward or reverse) match. For paired-end reads, both primary (not secondary alignments) read alignments must have matching 5 prime coordinates and share orientation. When a duplicate is detected the highest quality duplicate is kept and the others have the duplicate flag set.

![duplicate](media/duplicate.png)

The issue of when and whether to remove duplicates is debated. Duplicates as defined above may result from duplicate fragments made during PCR, which skew the results and should be removed. But duplicates may also be from two different fragments aligning to exact same coordinates, in which case their sequences may differ and they should be kept. Differentiating between these two types of duplicates is the problem.

Duplicates may be removed before alignment or after. The risk with removing duplicates before alignment is that you are actually removing real signal instead of duplicate fragments. The risk with removing duplicates after alignment is that you may remove fragments with unique sequences that happen to share coordinates with another alignment.

- `-r` Remove reads marked with duplicate flag.
- `--reference PATH` Path to reference genome used for alignment.

### samtools index

Indexes a coordinate sorted bgzipped compressed SAM, BAM, or CRAM file for random access.

### samtools stats

Provides summary statistics on BAM files and outputs them in text file. Many of these statistics come from values in the FLAG column. `samtools stats` counts the number of alignments for each samtools flag type and interprets the combination of certain flags with the goal of indicating how well the alignment went. Information on each of the flags is in the [Samtools specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

![flags](media/flag.png)

In the SAM or BAM file, the second column is FLAG. Each alignment has a flag value which is a unique combination of the flags in the table above. The number 77 for example is flags 1+4+8+64. `samtools stats` deduces the independent flags and creates file wide statistics on them.

When the mapping of a read is ambiguous, it may have multiple mappings. One mapping is considered primary, and all the others have the **secondary** flag set.

When the mapping of a read is chimeric, meaning non-linear, one piece is considered representative, and the other piece(s) are given the **supplementary** flag.

## Call Variants

### Strelka

- `--rna` (configure option) Applies specific settings for the rna-seq variant calling use case. This option is still in development.
- `--exome` (configure option) Applies settings for the whole exome seqquencing use case which include disabling high depth filters.
- `--callRegions PATH` (configure option) Pass in a bed file with regions to check for variants in. Bed file must be gzipped and tabixed.
- `--indelCandidates PATH` (configure option) Pass pre-predicted indels. If running the somatic workflow, its best practice to pass manta's output `candidateSmallIndels.vcf.gz` here. If running the germline workflow, this is not recommended.
- `--quiet` (run option) Sends error log to `${STRELKA_ANALYSIS_PATH}/workspace/pyflow.data/logs/pyflow_log.txt` instead of stdout.
- `--mode local` Runs locally as opposed to running on a cluster.

### Manta

- `--indelCandidates PATH`

## Annotate Variants

### SnpEff

Annotates variant with impact (high, moderate, low, or modifier), functional consequence (early stop codon, missense mutation, synonymous mutation, etc.), potential clinical significance, and more.

### bcftools concat

Concatenate two vcfs with the same sample set.

### bcftools annotate

- `rename-chrs PATH` Will rename the chromosomes in a vcf according to the text file passed.

### SnpSift

From the same author as SnpEff, SnpSift is a tool that annotates VCFs with specific information from other databases.

- `annotate` Will annotate VCF with rsids from another VCF.
