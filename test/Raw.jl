using Test: @test

using FASTQ

# ---- #

DA = joinpath(@__DIR__, "data", "Test")

# ---- #

for (di, le) in (("cDNA", 4), ("DNA", 4), ("MousecDNA", 6))

    di = joinpath(DA, di)

    fq_ = FASTQ.Raw.find(di)

    @test length(fq_) == le

    @test all([endswith(fq, ".gz") for fq in fq_])

    pa = joinpath(FASTQ.TE, "check_read")

    @info pa

    FASTQ.Raw.check_read(pa, fq_, 8)

    @test sum(collect(endswith(fi, "fastqc.html") for fi in readdir(pa))) == le

    @test isfile(pa, "multiqc_report.html")

    FASTQ.Raw.concatenate(fq_)

    @test length(readdir(joinpath(dirname(di), string(basename(di), "Concatenated")))) == 2

end

# ---- #

R1 = joinpath(DA, "DNA", "test_dna_4k.R1.fastq.gz")

R2 = replace(R1, "1" => "2")

SOR1 = joinpath(DA, "DNA", "test_dna_40k.R1.fastq.gz")

SOR2 = replace(SOR1, "1" => "2")

# ---- #

PA = joinpath(FASTQ.TE, "check_raw")

args = (FASTQ.TE, R1, R2, 8)

FASTQ.Raw.check_read(args...)

@test sum(collect(endswith(fi, "fastqc.html") for fi in readdir(PA))) == 2

FASTQ.Raw.check_read(args..., sor1 = SOR1, sor2 = SOR2)

@test sum(collect(endswith(fi, "fastqc.html") for fi in readdir(PA))) == 4

# ---- #

TR = joinpath(FASTQ.TE, "trim")

FASTQ.Raw.trim(TR, 8, R1, R2)

@test sum([endswith(fi, "fastq.gz") for fi in readdir(TR)]) == 2

@test sum([occursin("fastp", fi) for fi in readdir(TR)]) == 2

# ---- #

FASTQ.Raw.psuedoalign(ou, tr, n_jo, r1, r2, fr, sd)

# ---- #

FASTQ.Raw.align_cdna(pa, ge, n_jo, sa, r1, r2)

# ---- #

FASTQ.Raw.align_dna(pa, sa, ba, r1, r2, ge, n_jo, me)

# ---- #

FASTQ.Raw.align_single_cell_cdna()

# ---- #
