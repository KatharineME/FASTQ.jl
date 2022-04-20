module command

using JSON
using ..Fastq

include("apply_cdna_to_genome.jl")

include("apply_cdna_to_transcriptome.jl")

include("apply_germline_dna_to_genome.jl")

include("apply_somatic_dna_to_genome.jl")

include("concatenate_read.jl")

include("read_setting.jl")

end
