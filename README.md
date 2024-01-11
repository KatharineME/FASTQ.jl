Sequence analysis pipeline for raw DNA and cDNA reads :sparkles:

### [Use](#use)

### [Prepare Environment](#prepare-environment)

### [Get FASTQ.jl](#get-fastqjl)

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

# cDNA psuedoalignment and gene_x_sample creation
julia> FASTQ.Command.measure_gene_expression()
```

## Prepare Environment

#### 1. Run

Get [Homebrew](https://brew.sh).

```bash
brew install fastqc kallisto samtools bcftools

pip3 install multiqc
```

#### 2. Get SnpEff

[Download](http://pcingola.github.io/SnpEff/download/) and put in `tool/`.

#### 3. Get fastp, minimap2, and STAR

###### Mac

Unzip programs in `FASTQ.jl/tool/` and link their exectuables to `usr/local/bin/`

If programs in `tool/` fail then [compile programs](compile_program.md).

###### Linux

Get precompiled Linux versions from sources: [fastp](https://github.com/OpenGene/fastp), [minimap2](https://github.com/lh3/minimap2), [STAR](https://github.com/alexdobin/STAR).

#### 4. Get Docker, Manta, and Strelka

[Get docker](https://docs.docker.com/get-docker/).

Docker Desktop > Settings > Resources > File Sharing > Add `/tmp` and `/var`.

Download strelka-2.9.10.centos6_x86_64.tar.bz2 from [strelka releases](https://github.com/Illumina/strelka/releases).

Download manta-1.6.0.centos6_x86_64.tar.bz2 from [manta releases](https://github.com/Illumina/manta/releases).

Put manta and strelka in `tool/`.

#### 5. Set ulimit

Add `ulimit -n 10000` to `.zshrc`.

## Get FASTQ.jl

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

You can only p[ass the input directory which should only contain folders of samples. Each folder is the sample name. Inside are the reads for that sample.

All your tools are in one tool direcotry. This includes streak, manta, snpeff.

Your reference genome is in the same directory as all your chromosome files and your gtf / gif.

---

## :cowboy_hat_face: Howdy

To report a bug, request a feature, or leave a comment, just [submit an issue](https://github.com/GIT_USER_NAME/TEMPLATE.jl/issues/new/choose).

---

Powered by https://github.com/KwatMDPhD/PkgRepository.jl_
