using Test: @test, @test_throws

using FASTQ

# ---- #

const _RN1, _RN2 = FASTQ._RN1, FASTQ._RN2

const DA = joinpath(FASTQ._DA, "Test")

const DNA = joinpath(DA, "DNA")

const S1, S2 = "Sample1", "Sample2"

const FI1 = joinpath(DNA, S1, "test_dna_4k.R1.fastq.gz")

const FI2 = replace(FI1, _RN1 => _RN2)

const FI11 = joinpath(DNA, S2, "test_dna_40k.R1.fastq.gz")

const FI22 = replace(FI11, _RN1 => _RN2)

# ---- #

@test length(FASTQ.Support.check_os()) > 0

# ---- #

FASTQ.Support.start_docker()

@test length(readchomp(`docker ps -a`)) > 0

# ---- #

@test FASTQ.Support.calculate_size(FI1) > 286

# ---- #

@test FASTQ.Support.test_local_environment() === nothing

# ---- #

@test FASTQ.Support.log_top_level_function() === nothing

# ---- #

@test FASTQ.Support.log_sub_level_function() === nothing

# ---- #

@test FASTQ.Support.error_if_file_missing(readdir(joinpath(DNA, S1); join = true)) ===
      nothing

@test_throws ErrorException FASTQ.Support.error_if_file_missing((
    "Unicorn.txt",
    "Rainbow.txt",
))

# ---- #

const AB = joinpath("/Users", ENV["USER"], "Downloads")

# ---- #

@test FASTQ.Support.make_path_absolute("~/Downloads/") == AB

# ---- #

@test_throws ErrorException FASTQ.Support.make_sample_to_fastq_dictionary(DA, _RN1)

# ---- #

const TE = FASTQ.TE

# ---- #

TES1, TES2 = [FASTQ.Support.trash_remake_directory(joinpath(TE, pa)) for pa in (S1, S2)]

# ---- #

for (fi, sa) in ((FI1, TES1), (FI2, TES1), (FI11, TES2), (FI22, TES2))

    run(`cp $fi $(joinpath(sa, basename(fi)))`)

end

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
