module Fastq

include("run_docker_container.jl")

include("remove_docker_container.jl")

include("align_cdna.jl")

include("align_dna.jl")

include("call_germline_variant.jl")

include("call_somatic_variant.jl")

include("check_read.jl")

include("check_directory.jl")

include("concatenate.jl")

include("count_transcript.jl")

include("find.jl")

include("process_dna.jl")

include("process_soma_dna.jl")

include("test_environment.jl")

include("test_strelka_and_manta.jl")

include("trim.jl")

end
