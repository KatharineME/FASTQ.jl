include("find.jl")

n_jo, _, _, _, _, _, _, _, ou, _, _, _, _, _, _, _, _, _ = Fastq.read_setting(se)

Fastq.check_read(re_, joinpath(ou, "check_raw"), n_jo)
