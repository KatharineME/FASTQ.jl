module Fastq

include("align_cdna.jl")

include("align_dna.jl")

include("call_germline_variant.jl")

include("call_somatic_variant.jl")

include("check_read.jl")

include("concatenate.jl")

include("configure_and_run_manta.jl")

include("count_transcript.jl")

include("find.jl")

include("index_genome_files.jl")

include("make_directory.jl")

include("process_dna.jl")

include("process_soma_dna.jl")

include("reheader_vcf.jl")

include("remove_docker_container.jl")

include("run_docker_container.jl")

include("run_snpeff.jl")

include("test_local_environment.jl")

include("test_strelka_and_manta.jl")

include("trim.jl")

end
