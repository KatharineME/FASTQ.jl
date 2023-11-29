module Benchmark

using ..FASTQ

function _make_benchmark_files(vc, re, rt, nc, n_jo)

    red = splitext(re)[1]

    !isfile(red) ? run(`bgzip -d $re`) : nothing

    sd = replace(red, "fna" => "sdf")

    !isdir(sd) ? run(`$rt format -o $sd $red`) : nothing

    vqn = replace(vc, "pass" => "pass_rename_chromosomes")

    if !isfile(vqn)

        @warn "Renaming query vcf chromosomes"

        run(`bcftools annotate --threads=$n_jo --rename-chrs=$nc --output=$vqn $vc`)

    end

    red, sd, vqn

end


function _run_vcfeval(ou, vtr, vqn, sd, co, rt, n_jo)

    rte = joinpath(rt, "rtg")

    ouv = joinpath(ou, "Vcfeval")

    run(`$rte vcfeval 
        --baseline=$vtr 
        --bed-regions=$co 
        --calls=$vqn 
        --template=$sd 
        --output=$ouv
        --threads=$n_jo`)

    rte, ouv

end

function _run_happy(ha, vtr, vqn, sd, red, co, rt)

    ho, vo = "/home", "volume"

    vha = joinpath(ho, basename(ha[1]))

    pvt, pvqn, pbd, pre =
        [dirname(FASTQ.Support.make_path_absolute(pa)) for pa in (vtr, vqn, co, red)]

    vvt, vvqn, vre, vrt, vsd = [joinpath(ho, basename(pa)) for pa in (pvt, pvqn, pre, rt, sd)]

    vbd = joinpath(ho, "confident_regions_bed/")

    id = readchomp(`docker run 
            --interactive 
            --detach 
            --tty 
            --user root
            --$vo $pvt:$vvt 
            --$vo $pvqn:$vvqn 
            --$vo $pbd:$vbd 
            --$vo $pre:$vre 
            --$vo $ha:$vha 
            --$vo $rtg_tools:$vrt 
            --$vo $sd:$vsd 
            pkrusche/hap.py
            bash`)

    readchomp(
        `docker exec --interactive $id bash -c "/opt/hap.py/bin/hap.py $(joinpath(vvt, basename(vtr))) $(joinpath(vvqn, basename(vqn))) -f $(joinpath(vbd, basename(co))) -r $(joinpath(vre, basename(red))) -o $(joinpath(vha, "hap.py")) --engine-vcfeval-path $vrt --engine-vcfeval-template $vsd"`,
    )

end

function benchmark(
    output_directory,
    reference_genome,
    rtg_tools,
    number_of_jobs,
    query_vcf,
    truth_vcf,
    confident_regions_bed,
)
    FASTQ.Support.log_top_level_function()

    chs, chn = FASTQ.Reference.get_chromosome_file_path(reference_genome)
    FASTQ.Support.error_if_file_missing((chn, chs, query_vcf, truth_vcf, confident_regions_bed))

    ha = FASTQ.Support.make_analysis_directory(output_directory, "Benchmark", ("happy",))

    red, sd, vqn =
        _make_benchmark_files(query_vcf, reference_genome, rtg_tools, chn, number_of_jobs)

    @info "Running vcfeval"

    _run_vcfeval(
        output_directory,
        truth_vcf,
        vqn,
        sd,
        confident_regions_bed,
        rtg_tools,
        number_of_jobs,
    )

    FASTQ.Reference.index_genome_file(reference_genome, chs)

    @info "Running hap.py"

    _run_happy(ha, truth_vcf, vqn, sd, red, confident_regions_bed, rtg_tools)

    FASTQ.Support.remove_docker_container(id)

    nothing

end

end
