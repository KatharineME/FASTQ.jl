using Test: @test

using FASTQ

# ---- #

const TE = FASTQ.TE

const ID = mkdir(joinpath(TE, "Index"))

const GR = joinpath(FASTQ._DA, "GRCh38")

const GEN = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set"

const GEPA = joinpath(GR, GEN)

const GE = joinpath(GEPA, string(GEN, ".fna.gz"))

const CHS = joinpath(GEPA, "chromosome.bed.gz")

# ---- #

run(`cp $GE $CHS $ID`)

FASTQ.Reference.index_genome_file(joinpath(ID, basename(GE)), joinpath(ID, basename(CHS)))

@test lastindex([fi for fi in readdir(ID) if fi != ".DS_Store"]) == 5

# ---- #

const ST = mkdir(joinpath(TE, "GenerateStarFile"))

const GA = joinpath(GEPA, "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf")

# ---- #

run(`cp $GE $GA $ST`)

# ---- #

FASTQ.Reference.generate_star_genome_file(joinpath(ST, basename(GE)), 8; ga = nothing)

@test round(FASTQ.Support.calculate_size(joinpath(ST, "StarIndex", "SA"))) == 23
