include("_.jl")

n_jo, _, _, _, _, _, _, _, ou, ger1, ger2, _, _, _, _, _, _, _ = Fastq.read_setting(se)

Fastq.trim(ger1, ger2, joinpath(ou, "trim"), n_jo)
