module support

using ..Fastq

include("error_if_directory.jl")

include("get_full_path.jl")

include("index_genome_files.jl")

include("log.jl")

include("test_local_environment.jl")

include("test_strelka_and_manta.jl")

end
