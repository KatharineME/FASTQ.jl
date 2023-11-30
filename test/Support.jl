using Test: @test, @test_throws

using FASTQ

# ---- #

@test FASTQ.Support.test_local_environment() === nothing

# ---- #

@test FASTQ.Support.log_top_level_function() === nothing

# ---- #

@test FASTQ.Support.log_sub_level_function() === nothing

# ---- #

const DA = joinpath(FASTQ._DA, "Test")

const DNA = joinpath(DA, "DNA")

const S1 = "Sample1"

# ---- #

@test FASTQ.Support.error_if_file_missing(readdir(joinpath(DNA, S1), join = true)) === nothing

@test_throws ErrorException FASTQ.Support.error_if_file_missing(("Unicorn.txt", "Rainbow.txt"))

# ---- #

const AB = joinpath("/Users", ENV["USER"], "Downloads")

# ---- #

@test FASTQ.Support.make_path_absolute("~/Downloads/") == AB

# ---- #

const S2 = replace(S1, "1" => "2")

const _RN1 = FASTQ._RN1

# ---- #

@test_throws ErrorException FASTQ.Support.make_sample_to_fastq_dictionary(DA, _RN1)

# ---- #

const TE = FASTQ.TE

# ---- #

run(`cp -r $DNA"/" $TE`)

# ---- #

sa_fq_ = FASTQ.Support.make_sample_to_fastq_dictionary(TE, _RN1)

@test all(occursin("Sample", basename(ke)) for ke in keys(sa_fq_))

@test all(lastindex(get(sa_fq_, sa, "")) == 2 for sa in keys(sa_fq_))

# ---- #

const AN = joinpath(FASTQ.TE, "Analysis")

# ---- #

FASTQ.Support.make_analysis_directory(FASTQ.TE, "Analysis", ("One", "Two"); sa_fq_ = sa_fq_)

@test readdir(joinpath(AN, "One")) == [S1, S2]

# ---- #

@test FASTQ.Support.test_strelka_and_manta(joinpath(FASTQ.PR, "tool")) === nothing
