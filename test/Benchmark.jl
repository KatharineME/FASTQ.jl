using Test: @test

using FASTQ

# ---- #

const TE = FASTQ.TE

const DA = FASTQ._DA

const DAR = joinpath(DA, "GRCh38")

const GEN = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set"

const GEPA = joinpath(DAR, GEN)

const GE = joinpath(GEPA, string(GEN, ".fna.gz"))

const CHN = joinpath(GEPA, "chrn_n.tsv")

const RT = joinpath(FASTQ.PR, "tool", "rtg-tools-3.11")

const DAT = joinpath(DA, "Test")

const VQ = joinpath(DAT, "VCF", "Sample1", "concat.vcf.gz")

const BE = joinpath(DAT, "Benchmark")

const VT = joinpath(BE, "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")

const BED =
    joinpath(BE, "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed.gz")

const N_JO = 8

# ---- #

@test FASTQ.Benchmark.benchmark(TE, GE, RT, N_JO, VQ, VT, BED) === nothing

@test lastindex([
    fi for
    fi in readdir(joinpath(TE, "Benchmark", "Happy")) if fi != ".DS_Store"
]) == 11
