using Test: @test

using Nucleus

using FASTQ

# ---- #

# ----------------------------------------------------------------------------------------------- #

const FA = "FASTQ"

@test isconst(FASTQ, :PR)

basename(FASTQ.PR) == "$FA.jl"

# ---- #

@test isconst(FASTQ, :_DA)

@test basename(FASTQ._DA) == "data"

@test [fi for fi in readdir(FASTQ._DA) if fi != ".DS_Store"] == ["CellRangerBarcodes", "GRCh37", "GRCh38", "GRCm38", "Test"]

# ---- #

@test isconst(FASTQ, :TE)

@test basename(FASTQ.TE) == FA

@test isempty(readdir(FASTQ.TE))

# ---- #

const CE = "centos"

@test isconst(FASTQ, :_MA)

@test all(occursin(ne, FASTQ._MA) for ne in ("manta", CE))

# ---- #

@test isconst(FASTQ, :_ST)

@test all(occursin(ne, FASTQ._ST) for ne in ("strelka", CE))

# ---- #

const SR = joinpath(dirname(@__DIR__), "src")

const IG_ = (r"^[!_]",)

const MO_ = Nucleus.Path.read(SR; ig_ = IG_)

const TE_ = Nucleus.Path.read(@__DIR__; ig_ = IG_)

# ---- #

for jl in MO_

    @test chop(jl; tail = 3) == chop(readline(joinpath(SR, jl)); head = 7, tail = 0)

end

# ---- #

# @test symdiff(MO_, TE_) == ["FASTQ.jl", "runtests.jl"]

# ---- #

for jl in TE_

    if jl != "runtests.jl"

        @info "Testing $jl"

        run(`julia --project $jl`)

    end

end
