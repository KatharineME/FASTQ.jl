# DNA and cDNA sequence analysis

🚧 This README is under construction 🚧

## [Use](#use)

### [Prepare Environment](#prepare-environment)

### [Get FASTQ](#get-fastq)

### [Run Tests](#run-tests)

### [Programs Called](program.md)

## Use

```bash
cd FASTQ.jl

julia --project

julia> using FASTQ

# Concatenate fastq files
julia> FASTQ.command.concatenate_fastq()

# Germline DNA alignment and variant calling
julia> FASTQ.Command.call_variants_on_germline_dna()

# Somatic DNA alignment and variant calling
julia> FASTQ.Command.call_variants_on_somatic_dna()

# cDNA alignment and variant calling
julia> FASTQ.Command.call_variants_on_cdna()

# cDNA psuedoalignment
julia> FASTQ.Command.measure_gene_expression()
```

## Prepare Environment

#### 1. Run

Get [Homebrew](https://brew.sh).

```bash
brew install fastp fastqc kallisto samtools bcftools

pip3 install multiqc
```

#### 2. Get SnpEff

[Download](http://pcingola.github.io/SnpEff/download/) and put in `tool/`.

#### 3. Get minimap2 and STAR

###### Mac

Unzip programs in `FASTQ.jl/tool/` and link their exectuables to `usr/local/bin/`

If programs in `tool/` fail then [compile programs](compile_program.md).

###### Linux

Get precompiled Linux versions from sources: [minimap2](https://github.com/lh3/minimap2), [STAR](https://github.com/alexdobin/STAR).

#### 4. Get Docker, Manta, and Strelka

[Get docker](https://docs.docker.com/get-docker/).

Docker Desktop > Settings > Resources > File Sharing > Add `/tmp` and `/var`.

Download strelka-2.9.10.centos6_x86_64.tar.bz2 from [strelka releases](https://github.com/Illumina/strelka/releases).

Download manta-1.6.0.centos6_x86_64.tar.bz2 from [manta releases](https://github.com/Illumina/manta/releases).

Put manta and strelka in `tool/`.

#### 5. Set ulimit

Add `ulimit -n 10000` to `.zshrc`.

## Get FASTQ

```bash
git clone https://github.com/KatharineME/FASTQ.jl

cd FASTQ.jl

julia --project --eval "using Pkg; Pkg.instantiate()"

```

## Test

```bash
julia --project --eval "using Pkg; Pkg.test()"
```

## Keep in Mind

FASTQ.jl assumes the following

`input/` should be like this

```md
[ 160] input
├── [ 128] Sample1
│   ├── [ 0] R1.fastq.gz
│   └── [ 0] R2.fastq.gz
└── [ 128] Sample2
├── [ 0] R1.fastq.gz
└── [ 0] R2.fastq.gz
```

`tool/` should contain needed tools of correct versions

```md
[ 416] tool/
├── [ 544] STAR-2.7.9a
├── [ 15M] STAR-2.7.9a.zip
├── [ 192] manta-1.6.0.centos6_x86_64
├── [2.6K] minimap2-2.24
├── [1.1M] minimap2-2.24.zip
├── [ 416] rtg-tools-3.11
├── [5.1M] rtg-tools-3.11-nojre.zip
├── [ 352] snpEff
└── [ 192] strelka-2.9.10.centos6_x86_64
```

The reference genome should be adjacent to chromosome and gtf files.

```md
[ 576] GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set
├── [1.4G] GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf
├── [ 51M] GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
├── [2.9G] GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
├── [120K] GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai
├── [847M] GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
├── [120K] GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz.fai
├── [754K] GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz.gzi
├── [ 0] GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz.kallisto_index
├── [6.8G] GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz.mmi
├── [ 736] GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.sdf
├── [ 576] StarIndex
├── [ 202] chrn_n.tsv
├── [ 426] chromosome.bed
├── [ 252] chromosome.bed.gz
├── [4.1K] chromosome.bed.gz.tbi
└── [ 227] n_chrn.tsv
```

---

Made by [Kata](https://github.com/KwatMDPhD/Kata.jl) 🥋
