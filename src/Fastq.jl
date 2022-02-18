module Fastq

const MANTA = "manta-1.6.0.centos6_x86_64"

const STRELKA = "strelka-2.9.10.centos6_x86_64"

const TRIMMED_R1 = "trimmed.R1.fastq.gz"

const TRIMMED_R2 = "trimmed.R2.fastq.gz"

include("support/align_cdna.jl")

include("support/align_dna.jl")

include("support/call_germline_variant.jl")

include("support/call_somatic_variant.jl")

include("support/check_read.jl")

include("support/combine_vcf.jl")

include("support/concatenate.jl")

include("support/configure_and_run_manta.jl")

include("support/constant.jl")

include("support/find.jl")

include("support/get_full_path.jl")

include("support/index_genome_files.jl")

include("support/make_directory.jl")

include("support/psuedoalign.jl")

include("support/reheader_vcf.jl")

include("support/remove_docker_container.jl")

include("support/run_docker_container.jl")

include("support/run_snpeff.jl")

include("support/set_strelka_manta_run.jl")

include("support/test_local_environment.jl")

include("support/test_strelka_and_manta.jl")

include("support/trim.jl")

include("command/examine_read.jl")

include("command/process_dna.jl")

include("command/process_soma_dna.jl")

include("command/read_setting.jl")

end
