module vcf

using ..Fastq

include("annotate_with_snpeff.jl")

include("annotate_with_snpsift.jl")

include("combine_vcf.jl")

include("reheader_vcf.jl")

end
