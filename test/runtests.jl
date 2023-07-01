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
