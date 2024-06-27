module HLAType

using ..FASTQ

function hlatype(output_directory, bam, sample, tool_directory)

    FASTQ.Support.log_top_level_function()

    ba = basename(bam)

    bi = string(ba, ".bai")

    id = joinpath(dirname(bam), bi)

    FASTQ.Support.error_if_file_missing((bam, id))

    hl = joinpath(tool_directory, "HLA")

    ts = joinpath(hl, "tests")

    run(`cp -f $bam $ip $ts`)

    di = pwd()

    cd(hl)

    run(
        `docker run -v $hl:$hl -w $hl humanlongevity/hla --sample_id $sample --input_bam_path tests/$ba --output_path test`,
    )

    run(`cp -Rf $(joinpath(hl, string("hla-", sample))) $output_directory`)

    run(`cp -f $(joinpath(hl, "test", string("report-", sample, "-hla.json"))) $output_directory`)

    rm(joinpath(ts, ba))

    rm(joinpath(ts, bi))

    cd(di)

end

end
