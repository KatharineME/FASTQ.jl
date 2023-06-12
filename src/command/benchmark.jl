function benchmark(
    output_directory,
    reference_genome,
    rtg_tools,
    number_of_jobs,
    name_chromosome,
    query_vcf,
    truth_vcf,
    confident_regions_bed,
)

    Fastq.support.log()

    pa = joinpath(output_directory, "benchmark")

    Fastq.support.error_if_directory(pa)


    # Make vcfeval sdf

    red = split(reference_genome, ".gz")[1]

    if !isfile(reference_genome)

        run(`bgzip -d $reference_genome`)

    end

    sd = replace(red, "fna" => "sdf")

    if !isdir(sd)

        println("\nMaking vcfeval genome sdf\n")

        run(`$rtg_tools format -o $sd $red`)

    end


    # Rename chromosomes from numbers to strings

    vqn = replace(query_vcf, "pass" => "pass_rename_chromosomes")

    if !isfile(vqn)

        println("\nRenaming query VCF chromosomes\n")

        run(
            `bcftools annotate --threads=$number_of_jobs --rename-chrs=$name_chromosome --output=$vqn $query_vcf`,
        )

    end


    # Run vcfeval

    println("\nRunning vcfeval\n")

    ouv = joinpath(pa, "vcfeval")

    rte = joinpath(rtg_tools, "rtg")

    run(`$rte vcfeval 
        --baseline=$truth_vcf 
        --bed-regions=$confident_regions_bed 
        --calls=$vqn 
        --template=$sd 
        --output=$ouv
        --threads=$number_of_jobs`)


    # Run hap.py in container

    println("\nRunning hap.py\n")

    ouh = joinpath(pa, "happy/")

    mkdir(ouh)

    ho = "/home/"

    vouh = joinpath(ho, splitpath(ouh)[end])

    pvt = dirname(Fastq.support.get_full_path(truth_vcf))

    vvt = joinpath(ho, basename(pvt))

    pvqn = dirname(Fastq.support.get_full_path(vqn))

    vvqn = joinpath(ho, basename(pvqn))

    pbd = dirname(Fastq.support.get_full_path(confident_regions_bed))

    vbd = joinpath(ho, "confident_regions_bed/")

    pre = dirname(Fastq.support.get_full_path(red))

    vre = joinpath(ho, basename(pre))

    vrt = joinpath(ho, basename(rtg_tools))

    vsd = joinpath(ho, basename(sd))

    id = readlines(pipeline(`docker run 
            --interactive 
            --detach 
            --tty 
            --user root
            -v $pvt:$vvt 
            -v $pvqn:$vvqn 
            -v $pbd:$vbd 
            -v $pre:$vre 
            -v $ouh:$vouh 
            -v $rtg_tools:$vrt 
            -v $sd:$vsd 
            pkrusche/hap.py
            bash`))

    readlines(
        pipeline(
            `docker exec --interactive $id bash -c "/opt/hap.py/bin/hap.py $(joinpath(vvt, basename(truth_vcf))) $(joinpath(vvqn, basename(vqn))) -f $(joinpath(vbd, basename(confident_regions_bed))) -r $(joinpath(vre, basename(red))) -o $(joinpath(vouh, "hap.py")) --engine-vcfeval-path $vrt --engine-vcfeval-template $vsd"`,
        ),
    )

    Fastq.support.remove_docker_container(id)

end
