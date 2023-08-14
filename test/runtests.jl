using Logging

using Test

using BioLab: BioLab, @is_error

using FASTQ

TE = FASTQ.TE














using Aqua

include("environment.jl")

# ---- #

# Aqua.test_all(FASTQ; ambiguities = false)
# 
# Aqua.test_ambiguities(FASTQ)

# ----------------------------------------------------------------------------------------------- #

te_ = filter!(!startswith('_'), readdir(@__DIR__))

# ---- #

for te in te_

    if te != "runtests.jl"

        @info "Testing $te"

        run(`julia --project $te`)

    end

end


















using Test: @test

using BioLab

# ---- #

# ----------------------------------------------------------------------------------------------- #

@test isconst(BioLab, :_DA)

@test basename(BioLab._DA) == "data"

@test BioLab.Path.read(BioLab._DA) == [
    "CLS",
    "DataFrame",
    "Dict",
    "FeatureSetEnrichment",
    "GCT",
    "GMT",
    "Gene",
    "Plot",
    "SingleCell",
]

# ---- #

@test isconst(BioLab, :TE)

@test basename(BioLab.TE) == "BioLab"

@test isempty(BioLab.Path.read(BioLab.TE))

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

@test symdiff(MO_, TE_) == ["BioLab.jl", "runtests.jl"]

# ---- #

for jl in TE_

    if jl != "runtests.jl"

        @info "Testing $jl"

        run(`julia --project $jl`)

    end

endd
