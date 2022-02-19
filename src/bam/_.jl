module bam

include("call_germline_variant.jl")

include("call_somatic_variant.jl")

include("configure_and_run_manta.jl")

include("remove_docker_container.jl")

include("run_docker_container.jl")

include("run_snpeff.jl")

include("set_strelka_manta_run.jl")

end
