#=
Illumina TruSeq adapters:

--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA 
--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
=#

function trim(
    fq1::String,
    fq2::String,
    pa::String,
    n_jo::Int,
)::Nothing

    println("Trimming...")

    mkpath(pa)
    
    ou1 = joinpath(pa, basename(fq1))

    ou2 = joinpath(pa, basename(fq2))

    run(`fastp --detect_adapter_for_pe --html fastp.html --in1 $fq1 --in2 $fq2 --out1 $ou1 --out2 $ou2`)

    return nothing

end
