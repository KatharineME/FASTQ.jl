module FASTQ

const ROOT = pkgdir(Fastq)

const MANTA = "manta-1.6.0.centos6_x86_64"

const STRELKA = "strelka-2.9.10.centos6_x86_64"

const TRIMMED_R1 = "trimmed.R1.fastq.gz"

const TRIMMED_R2 = "trimmed.R2.fastq.gz"

include("Abundance.jl")

include("BAM.jl")

include("Command.jl")

include("Raw.jl")

include("Support.jl")

include("VCF.jl")

end
