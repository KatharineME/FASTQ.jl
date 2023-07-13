Sequence analysis pipeline for raw DNA and cDNA reads :sparkles:

### [Use](#use)

### [Prepare Environment](#prepare-environment)

### [Get Fastq.jl](#get-fastqjl)

### [Run Tests](#run-tests)

### [Set settings.json](setting.md)

### [Programs Called](program.md)

## Use

Copy `Fastq.jl/setting.json` to your project directory and edit it accordingly. See [setting.md](setting.md) for details.

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

Update `Fastq.jl/test/setting.json`.

```bash
julia --project --eval "using Pkg; Pkg.test()"
```

---

## :cowboy_hat_face: Howdy

To report a bug, request a feature, or leave a comment, just [submit an issue](https://github.com/GIT_USER_NAME/TEMPLATE.jl/issues/new/choose).

---

Powered by https://github.com/KwatMDPhD/PkgRepository.jl_
