module Fastq

include("support/align_cdna.jl")

include("support/align_dna.jl")

include("support/call_germline_variant.jl")

include("support/call_somatic_variant.jl")

include("support/check_read.jl")

include("support/combine_vcf.jl")

include("support/concatenate.jl")

include("support/configure_and_run_manta.jl")

include("support/constant.jl")

include("support/count_transcript.jl")

include("support/find.jl")

include("support/get_full_path.jl")

include("support/index_genome_files.jl")

include("support/make_directory.jl")

include("support/reheader_vcf.jl")

include("support/remove_docker_container.jl")

include("support/run_docker_container.jl")

include("support/run_snpeff.jl")

include("support/test_local_environment.jl")

include("support/test_strelka_and_manta.jl")

include("support/trim.jl")

include("workflow/examine_read.jl")

include("workflow/process_dna.jl")

include("workflow/process_soma_dna.jl")

include("workflow/read_setting.jl")

end
