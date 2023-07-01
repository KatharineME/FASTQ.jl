module FASTQ

const DA = joinpath(dirname(@__DIR__), "data")

const TE = joinpath(tempdir(), "FASTQ")

const MA = "manta-1.6.0.centos6_x86_64"

const ST = "strelka-2.9.10.centos6_x86_64"

const TR1 = "trimmed.R1.fastq.gz"

const TR2 = "trimmed.R2.fastq.gz"

for jl in readdir(@__DIR__)

    if !startswith(jl, '_') && jl != "FASTQ.jl"

        include(jl)

    end

end

macro is_error(ex)

    quote

        try

            $(esc(ex))

            false

        catch er

            @info "Errored." er

            true

        end

    end

end

function __init__()

    rm(TE; recursive = true, force = true)

    mkdir(TE)

end

end
