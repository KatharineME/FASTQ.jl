function benchmark(se)

    Fastq.support.log()

    fe_va = read_setting(se)

    ou, re, rt, n_jo, ch, vq, vt, bd = fe_va["output_directory"],
    fe_va["reference_genome"],
    fe_va["rtg_tools"],
    fe_va["number_of_jobs"],
    fe_va["name_chromosome"],
    fe_va["query_vcf"],
    fe_va["truth_vcf"],
    fe_va["confident_regions_bed"]

    pa = joinpath(ou, "benchmark")

    Fastq.support.error_if_directory(pa)


    # Make vcfeval sdf

    red = split(re, ".gz")[1]

    if !isfile(re)

        run(`bgzip -d $re`)

    end

    sd = replace(red, "fna" => "sdf")

    if !isdir(sd)

        println("\nMaking vcfeval genome sdf\n")

        run(`$rt format -o $sd $red`)

    end


    # Rename chromosomes from numbers to strings

    vqn = replace(vq, "pass" => "pass_rename_chromosomes")

    if !isfile(vqn)

        println("\nRenaming query VCF chromosomes\n")

        run(`bcftools annotate --threads=$n_jo --rename-chrs=$ch --output=$vqn $vq`)

    end


    # Run vcfeval

    println("\nRunning vcfeval\n")

    ouv = joinpath(pa, "vcfeval")

    rte = joinpath(rt, "rtg")

    run(`$rte vcfeval 
        --baseline=$vt 
        --bed-regions=$bd 
        --calls=$vqn 
        --template=$sd 
        --output=$ouv
        --threads=$n_jo`)


    # Run hap.py in container

    println("\nRunning hap.py\n")

    ouh = joinpath(pa, "happy/")

    mkdir(ouh)

    ho = "/home/"

    vouh = joinpath(ho, splitpath(ouh)[end])

    pvt = dirname(Fastq.support.get_full_path(vt))

    vvt = joinpath(ho, basename(pvt))

    pvqn = dirname(Fastq.support.get_full_path(vqn))

    vvqn = joinpath(ho, basename(pvqn))

    pbd = dirname(Fastq.support.get_full_path(bd))

    vbd = joinpath(ho, "confident_regions_bed/")

    pre = dirname(Fastq.support.get_full_path(red))

    vre = joinpath(ho, basename(pre))

    vrt = joinpath(ho, basename(rt))

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
            -v $rt:$vrt 
            -v $sd:$vsd 
            pkrusche/hap.py
            bash`))

    readlines(
        pipeline(
            `docker exec --interactive $id bash -c "/opt/hap.py/bin/hap.py $(joinpath(vvt, basename(vt))) $(joinpath(vvqn, basename(vqn))) -f $(joinpath(vbd, basename(bd))) -r $(joinpath(vre, basename(red))) -o $(joinpath(vouh, "hap.py")) --engine-vcfeval-path $vrt --engine-vcfeval-template $vsd"`,
        ),
    )

    Fastq.support.remove_docker_container(id)

end
