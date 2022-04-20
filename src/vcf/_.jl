module vcf

using ..Fastq

include("combine_vcf.jl")

include("reheader_vcf.jl")

include("run_snpeff.jl")

end
