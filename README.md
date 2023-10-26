Sequence analysis pipeline for raw DNA and cDNA reads :sparkles:

### [Use](#use)

### [Prepare Environment](#prepare-environment)

### [Get Fastq.jl](#get-fastqjl)

### [Run Tests](#run-tests)

### [Set settings.json](setting.md)

### [Programs Called](program.md)

## Use

```bash
cd Fastq.jl

julia --project

julia> using Fastq

# Concatenate fastq files
julia> Fastq.command.concatenate_fastq("my/project/setting.json")

# Germline DNA alignment and variant calling
julia> FASTQ.Command.call_variants_on_germline_dna("my/project/setting.json")

# Somatic DNA alignment and variant calling
julia> FASTQ.Command.call_variants_on_somatic_dna("my/project/setting.json")

# cDNA alignment and variant calling
julia> FASTQ.Command.call_variants_on_cdna("my/project/setting.json")

# cDNA psuedoalignment and gene_x_sample creation
julia> FASTQ.Command.measure_gene_expression("my/project/setting.json")
```

## Prepare Environment

_Currently works only on MacOS._

#### 1. Run

```bash
brew install fastqc kallisto samtools bcftools

pip3 install multiqc
```

#### 2. Download SnpEff

Download from [here](http://pcingola.github.io/SnpEff/download/) and link to `usr/local/bin`.

#### 3. Unzip programs in `Fastq.jl/tool/` and link their exectuables to `usr/local/bin/`

If programs in `tool/` fail then [compile programs](compile_program.md).

#### 4. Get Docker, Manta, and Strelka

[Get docker](https://docs.docker.com/get-docker/).

Download strelka-2.9.10.centos6_x86_64.tar.bz2 from [strelka releases](https://github.com/Illumina/strelka/releases).

Download manta-1.6.0.centos6_x86_64.tar.bz2 from [manta releases](https://github.com/Illumina/manta/releases).

Put manta and strelka in the same directory and unzip each. This directory is the "tool_directory" in `setting.json`.

#### 5. Set ulimit

Add `ulimit -n 10000` to `.zshrc` or `.zprofile`.

## Get Fastq.jl

```bash
git clone https://github.com/KatharineME/Fastq.jl

cd Fastq.jl

julia --project --eval "using Pkg; Pkg.instantiate()"

```

## Test

```bash
julia --project --eval "using Pkg; Pkg.test()"
```

## Old settings.json

{
"number_of_jobs": 8,
"memory": 8,
"organism": "human",
"molecule": "dna",
"exome": false,
"fragment_length": 51,
"fragment_length_standard_deviation": 0.05,
"sample": "test",
"annotate_with_rsid": true,
"tool_directory": "/Users/kate/craft/jl/Fastq.jl/tool",
"output_directory": "/Users/kate/craft/jl/Fastq.jl/test/data/output",
"read_name_scheme": "R1",
"read1": "/Users/kate/craft/jl/Fastq.jl/test/data/dna/test_dna_40k.R1.fastq.gz",
"read2": "/Users/kate/craft/jl/Fastq.jl/test/data/dna/test_dna_40k.R2.fastq.gz",
"somatic_read1": "/Users/kate/craft/jl/Fastq.jl/test/data/dna/test_dna_4k.R1.fastq.gz",
"somatic_read2": "/Users/kate/craft/jl/Fastq.jl/test/data/dna/test_dna_4k.R2.fastq.gz",
"dna_read_directory": "/Users/kate/craft/jl/Fastq.jl/test/data/dna/",
"cdna_read_directory": "/Users/kate/craft/jl/Fastq.jl/test/data/cdna/",
"reference_genome": "/Users/kate/craft/data/grch/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz",
"reference_transcriptome": "/Users/kate/craft/data/mouse_reference_transcriptome/Mus_musculus.GRCm38.cdna.all.fa.gz",
"mouse_transcript_to_mouse_gene": "/Users/kate/craft/data/mouse_transcript_mouse_gene.tsv",
"chromosome_position": "/Users/kate/craft/data/grch/chromosome/chromosome.bed.gz",
"chromosome_name": "/Users/kate/craft/data/grch/chromosome/chrn_n.tsv",
"name_chromosome": "/Users/kate/craft/data/grch/chromosome/n_chrn.tsv",
"variant_database": "/Users/kate/craft/Fastq.jl/tool/ensembl/homo_sapiens-chr1_y.vcf.gz",
"snpeff": "/Users/kate/craft/jl/Fastq.jl/tool/snpEff/snpEff.jar",
"truth_vcf": "/Users/kate/craft/data/benchmark/HG002_truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
"query_vcf": "/Users/kate/craft/data/benchmark/apply_germline_dna_to_genome/call_germline_variant/pass.vcf.gz",
"confident_regions_bed": "/Users/kate/craft/data/benchmark/HG002_truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed.gz",
"rtg_tools": "/Users/kate/craft/jl/Fastq.jl/tool/rtg-tools-3.11"
}

## Old settings.json Instructions

`number_of_jobs` _Int_<br>
Number of jobs

`memory` _Int_<br>
Gigabytes of memory

`organism` _String_<br>
"human" or "mouse" used by `apply_cdna_to_transcriptome` when creating gene_x_sample matrix

`molecule` _String_<br>
"dna" or "cdna"

`exome` _Bool_<br>
"true" or "false"

`fragment_length` _Int_<br>
The most common cDNA fragment length in your data, used by `apply_cdna_to_genome`

`fragment_length_standard_deviation` _Int_<br>
The estimated standard deviation of fragment lengths, used by `apply_cdna_to_genome`

`sample` _String_<br>
Sample name used for output file names and name of sample column in VCF

`annotate_with_rsid` _Bool_<br>
"true" or "false"

`tool_directory` _String_<br>
Path to directory storing strelka and manta used by `apply_germline_dna_to_genome`, `apply_somatic_dna_to_genome`, and `apply_cdna_to_genome`

`output_directory` _String_<br>
Path to where you want results outputted

`read_name_scheme` _String_<br>
Naming scheme your reads use, "R1" or "read1" for example

`read1` _String_<br>
Path to DNA forward read file used by `apply_germline_dna_to_genome` and `apply_somatic_dna_to_genome`

`read2` _String_<br>
Path to DNA reverse read file used by `apply_germline_dna_to_genome` and `apply_somatic_dna_to_genome`

`somatic_read1` _String_<br>
Path to somatic DNA forward read file used by `apply_somatic_dna_to_genome`

`somatic_read2` _String_<br>
Path to somatic DNA reverse read file used by `apply_somatic_dna_to_genome`

`dna_read_directory` _String_<br>
Path to directory with DNA read files used by `concatenate_fastq` to combine fastqs of the same read direction

`cdna_read_directory` _String_<br>
Path to directory with cDNA read files used by `apply_cdna_to_genome` and `apply_cdna_to_transcriptome` which both respect read file name and directory structure to keep different cDNA samples separate

`reference_genome` _String_<br>
Path to human reference genome fna.gz or fa.gz file used by `apply_germline_dna_to_genome`, `apply_somatic_dna_to_genome`, and `apply_cdna_to_genome`

`reference_transcriptome` _String_<br>
Path to human reference transcriptome fna.gz or fa.gz file used by `apply_cdna_to_transcriptome`

`mouse_transcript_to_mouse_gene` _String_<br>
Path to tsv with mouse transcript name column and mouse gene name column used by `apply_cdna_to_transcriptome` when `organism` is set to "mouse"

`chromosome_position` _String_<br>
Path to tsv with 3 columns: chromosome name, start position of chromosome (0), and end position of chromosome. Used by `apply_germline_dna_to_genome`, `apply_somatic_dna_to_genome`, and `apply_cdna_to_genome`

`chromosome_name` _String_<br>
Path to tsv with column of chromosome names ("chr1", "chr2") followed by a column of their integer names (1, 2)

`name_chromosome` _String_<br>
Path to tsv with column of chromosome integer names (1, 2) followed by a column of their names ("chr1", "chr2"), used by `benchmark`

`variant_database` _String_<br>
Path to variant database VCF file (dbsnp or ensembl variant database for example) used by `apply_germline_dna_to_genome`, `apply_somatic_dna_to_genome`, and `apply_cdna_to_genome` when `annotate_with_rsid` is set to "true"

`snpeff` _String_<br>
Path to snpEff.jar file used by `apply_germline_dna_to_genome`, `apply_somatic_dna_to_genome`, and `apply_cdna_to_genome`

`truth_vcf` _String_<br>
Path to truth VCF file (HG002_GRCh38_1_22_v4.2.1_benchmark.vcf for example) used by `benchmark`

`query_vcf`
Path to query VCF file containing variant calls generated by your workflow that you want to benchmark, used by `benchmark`

`confident_regions_bed`
Bed file listing highly confident regions where if a variant within these regions is missed it is called a false negative, used by `benchmark`

`rtg_tools`
Path to rtg-tools directory containing rtg executable, used by `benchmark`

---

## :cowboy_hat_face: Howdy

To report a bug, request a feature, or leave a comment, just [submit an issue](https://github.com/GIT_USER_NAME/TEMPLATE.jl/issues/new/choose).

---

Powered by https://github.com/KwatMDPhD/PkgRepository.jl_
