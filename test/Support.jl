using Test: @test, @test_throws

using FASTQ

# ---- #

@test FASTQ.Support.test_local_environment() == nothing

# ---- #

@test FASTQ.Support.log_top_level_function() == nothing

# ---- #

@test FASTQ.Support.log_sub_level_function() == nothing

# ---- #

const DA = joinpath(FASTQ._DA, "Test")

const DNA = joinpath(DA, "DNA")

@test FASTQ.Support.error_if_file_missing(readdir(DNA, join = true)) == nothing

@test_throws ErrorException FASTQ.Support.error_if_file_missing(("Unicorn.txt", "Rainbow.txt"))

# ---- #

const AB = joinpath("/Users", ENV["USER"], "Downloads")

@test FASTQ.Support.make_path_absolute("~/Downloads/") == AB

# ---- #

const S1 = "Sample1"

const S2 = "Sample2"

@test_throws ErrorException FASTQ.Support.make_sample_to_fastq_dictionary(DA)

SA1, SA2, = [mkdir(joinpath(FASTQ.TE, sa)) for sa in (S1, S2)]

const CP_ = ((SA1, "40k"), (SA2, "4k"))

for fi in readdir(DNA, join = true)

    for (sa, st) in CP_

        if occursin(st, fi)

            cp(fi, joinpath(sa, basename(fi)))

        end

    end

end

const SA_FQ_ = FASTQ.Support.make_sample_to_fastq_dictionary(FASTQ.TE)

@test [keys(SA_FQ_)...] == [SA1, SA2]

for sa in (SA1, SA2)

    @test length(get(SA_FQ_, sa, "")) == 2

end

# ---- #

const AN = joinpath(FASTQ.TE, "Analysis")

FASTQ.Support.make_analysis_directory(FASTQ.TE, "Analysis", ("One", "Two"); sa_fq_ = SA_FQ_)

@test readdir(joinpath(AN, "One")) == [S1, S2]

# ---- #

FASTQ.Support.test_strelka_and_manta(joinpath(FASTQ.PR, "tool"))
