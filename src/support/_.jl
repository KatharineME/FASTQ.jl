module support

using ..Fastq

include("get_full_path.jl")

include("index_genome_files.jl")

include("error_if_directory.jl")

include("test_local_environment.jl")

include("test_strelka_and_manta.jl")

end
