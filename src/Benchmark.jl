module Benchmark

using ..FASTQ

function _make_benchmark_files(vc, re, rt, nc, n_jo)

    red = splitext(re)[1]

    if !isfile(red)

        run(`bgzip -d $re`)

    end

    sd = replace(red, "fna" => "sdf")

    if !isdir(sd)

        @info "Making vcfeval genome sdf"

        run(`$rt format -o $sd $red`)

    end

    vqn = replace(vc, "pass" => "pass_rename_chromosomes")

    if !isfile(vqn)

        @warn "Renaming query vcf chromosomes"

        run(`bcftools annotate --threads=$n_jo --rename-chrs=$nc --output=$vqn $vc`)

    end

    red, sd, vqn

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

    rte = joinpath(rtg_tools, "rtg")

    ouv = joinpath(output_directory, "Vcfeval")

    run(`$rte vcfeval 
        --baseline=$truth_vcf 
        --bed-regions=$confident_regions_bed 
        --calls=$vqn 
        --template=$sd 
        --output=$ouv
        --threads=$number_of_jobs`)

    FASTQ.Reference.index_genome_file(reference_genome, chs)

    @info "Running hap.py"

    ho = "/home"

    vha = joinpath(ho, basename(ha[1]))

    pvt, pvqn, pbd, pre = [
        dirname(FASTQ.Support.make_path_absolute(pa)) for
        pa in (truth_vcf, vqn, confident_regions_bed, red)
    ]

    vvt, vvqn, vre, vrt, vsd =
        [joinpath(ho, basename(pa)) for pa in (pvt, pvqn, pre, rtg_tools, sd)]

    vbd = joinpath(ho, "confident_regions_bed/")

    vo = "volume"

    id = readlines(pipeline(`docker run 
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
            bash`))


    vh = joinpath(vha, "hap.py")

    vr = joinpath(vre, basename(red))

    vb = joinpath(vbd, basename(confident_regions_bed))

    vtr = joinpath(vvt, basename(truth_vcf))

    vqu = joinpath(vvqn, basename(vqn))

    readlines(
        pipeline(
            `docker exec --interactive $id bash -c "/opt/hap.py/bin/hap.py $vtr $vqu -f $vb -r $vr -o $vh --engine-vcfeval-path $vrt --engine-vcfeval-template $vsd"`,
        ),
    )

    FASTQ.Support.remove_docker_container(id)

    nothing

end

end
