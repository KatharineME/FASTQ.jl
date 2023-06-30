module FASTQ

const DA = joinpath(dirname(@__DIR__), "data")

const TE = joinpath(tempdir(), "FASTQ")

const MA = "manta-1.6.0.centos6_x86_64"

const ST = "strelka-2.9.10.centos6_x86_64"

const TR1 = "trimmed.R1.fastq.gz"

const TR2 = "trimmed.R2.fastq.gz"

include("Abundance.jl")

include("BAM.jl")

include("Command.jl")

include("Raw.jl")

include("Support.jl")

include("VCF.jl")

end
