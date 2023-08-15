using Test: @test

using BioLab

using FASTQ

# ---- #

# ----------------------------------------------------------------------------------------------- #

# TODO: add constants
@test isconst(FASTQ, :_DA)

@test basename(FASTQ._DA) == "data"

@test BioLab.Path.read(FASTQ._DA) ==
      ["Clinvar", "Ensembl", "Gene", "Mouse", "ReferenceGenome", "ReferenceTranscriptome", "Test"]

# ---- #

@test isconst(FASTQ, :TE)

@test basename(FASTQ.TE) == "FASTQ"

@test isempty(BioLab.Path.read(FASTQ.TE))

# ---- #

const SR = joinpath(dirname(@__DIR__), "src")

const IG_ = (r"^[!_]",)

const MO_ = BioLab.Path.read(SR; ig_ = IG_)

const TE_ = BioLab.Path.read(@__DIR__; ig_ = IG_)

# ---- #

for jl in MO_

    @test chop(jl; tail = 3) == chop(readline(joinpath(SR, jl)); head = 7, tail = 0)

end

# ---- #

@test symdiff(MO_, TE_) == ["FASTQ.jl", "runtests.jl"]

# ---- #

for jl in TE_

    if jl != "runtests.jl"

        @info "Testing $jl"

        run(`julia --project $jl`)

    end

end
