module bam

using ..Fastq

include("call_germline_variant.jl")

include("call_somatic_variant.jl")

include("configure_and_run_manta.jl")

include("run_strelka_manta_docker_container.jl")

include("set_strelka_manta_run.jl")

end
