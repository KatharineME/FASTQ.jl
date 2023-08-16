using Test: @test

using FASTQ

# ---- #

DA = joinpath(@__DIR__, "data", "Test")

# ---- #

for (di, le) in (("cDNA", 4), ("DNA", 4), ("DNAConcatenated", 2), ("MousecDNA", 6))

    di = joinpath(DA, di)

    fq_ = FASTQ.Raw.find(di)

    @test length(fq_) == le

    @test all([endswith(fq, ".gz") for fq in fq_])

end

# ---- #


# ---- #
