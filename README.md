# Fastq.jl

Turns raw sequencing reads into intpretable files like: VCF for variants or a gene by sample matrix for gene expression analysis.

For each of the use cases below, there is a corresponding Jupyter Notebook workflow:

Germline DNA > VCF (Genome alignment)

Somatic DNA > VCF (Genome alignment)

cDNA > VCF (Genome alignment)

cDNA > Gene by sample (Psuedoalignment)

## Alignment

### minimap2

#### `-x ` 
A recommended meta flag that specifies alginment chain bandwidth, elongation, discard, how many secondary alignment should be output, and more

#### `-a` 
Generates CIGAR and output alignments in SAM format

#### `--sr` 
Enable short read mode

#### `--splice` 
Enable splice alignment mode

#### `-uf`
Use transcript strand to find canonical splice sites. 

#### `-t NUM`
Number of threads

#### `-K NUM`
Number of bases loaded into memory to process in a mini-batch [500M]. A large NUM pased here helps with load balancing in the multi-threading mode.

#### `-R`
SAM read group line in a format like @RG\\tID:foo\\tSM:bar [].

### samtools fixmate

Corrects flaws in read-pairing that the aligner may have introduced. It ensures the SAM FLAG, RNEXT, PNEXT, and TLEN fields are correct and consistent.

#### `-u`
Uncompressed output

#### `-m`
Add mate score tag

#### `-` 
A synonym for stdin and stdout in samtools

### samtools sort

Orders aligned reads by chromosome and coordinate. sort is highly parallel so adding threads is beneficial.

#### `-u`
Uncompressed output

#### `-l 1`
Fastest level of BAM compressed output

#### `--threads NUM`
Threads

#### `-T PATH`
Specificy a temporary storage directory like `/tmp/example_prefix/`

### samtools markdup

Marks duplicate alignments in a coordinate sorted file using flags added by `fixmate -m`. Duplicates are defined as aligned reads whose 5 prime coordinates and orientation (forward or reverse) match. For paried-end reads, both primary (not secondary alignments) read alignments must have matching 5 prime coordinates and share orientation. When a duplicate is detected the highest quality duplicate is kept and the others have the duplicate flag set.

The issue of when and whether to remove duplicates is debated. Duplicates as defined above may result from duplicate fragments made during PCR, which skew the results and should be removed. But duplicates may also be from two different fragments aligning to exact same corrdinates, in which case their sequences may differ and they should be kept. Differentiating between these two types of duplicates is the problem.  

Duplicates may be removed before alignment or after. The risk with removing duplicates before alignment is that you are actually removing real signal instead of duplicate fragments. The risk with removing duplicates after alignment is that you may remove fragments with unique sequences that happen to share coordinates with another alignment.

#### `-r`
Remove reads marked with duplicate flag

#### `-@NUM`
Threads

#### `--reference PATH`
Path to reference genome used for alignment

### samtools index

Indexes a coordinate sorted bgzipped compressed SAM, BAM, or CRAM file for random access.

#### `--threads NUM`
Threads

### samtools flagstat

Counts the number of alignments for each samtools flag type and indicates how successful the alignment was. Information on each of the flags is in the [Samtools specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

When the mapping of a read is ambiguous, it may have multiple mappings. One mapping is considered primary, and all the others have the __secondary__ flag set.

When the mapping of a read is chimeric, meaning non-linear, one piece is considered representative, and the other piece(s) are given the __supplementary__ flag.

#### `--threads NUM`
Threads

#### `--output-fmt FORMAT`
Sets output format to "tsv" or "json"

## Contribution

To report a bug, request a feature, or leave a comment (about anything related to this repository), just [submit an issue](https://github.com/KatharineME/Fastq.jl.jl/issues/new/choose).

---

Made by https://github.com/KwatMDPhD/PkgRepository.jl
