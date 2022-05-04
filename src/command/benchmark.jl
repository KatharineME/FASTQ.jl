function benchmark(se)

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


    # Make vcfeval SDF

    red = split(re, ".gz")[1]

    if !isfile(re)

        run(`bgzip -d $re`)

    end

    sd = replace(red, "fna" => "sdf")

    if !isdir(sd)

        run(`$rt format -o $sd $red`)

    end


    # Rename chromosomes from numbers to strings

    vqn = replace(vq, "pass" => "pass_rename_chromosomes")

    if !isfile(vqn)

        run(`bcftools annotate --threads=$n_jo --rename-chrs=$ch --output=$vqn $vq`)

    end


    # run vcfeval

    ouv = joinpath(pa, "vcfeval")

    rte = joinpath(rt, "rtg")

    run(`$rte vcfeval 
        --baseline=$vt 
        --bed-regions=$bd 
        --calls=$vqn 
        --template=$sd 
        --output=$ouv
        --threads=$n_jo`)


    # run hap.py in hap.py docker container

    ouh = joinpath(pa, "happy/")

    mkdir(ouh)

    ho = "/home/"

    vouh = joinpath(ho, splitpath(ouh)[end])

    println(vouh)

    pvt = dirname(Fastq.support.get_full_path(vt))

    vvt = joinpath(ho, basename(pvt))

    println(pvt)

    println(vvt)

    pvqn = dirname(Fastq.support.get_full_path(vqn))

    vvqn = joinpath(ho, basename(pvqn))

    println(pvqn)

    println(vvqn)

    pbd = dirname(Fastq.support.get_full_path(bd))

    vbd = joinpath(ho, "confident_regions_bed/")

    println(pbd)

    println(vbd)

    pre = dirname(Fastq.support.get_full_path(red))

    vre = joinpath(ho, basename(pre))

    println(pre)

    println(vre)

    vrt = joinpath(ho, basename(rt))

    println(vrt)

    vsd = joinpath(ho, basename(sd))

    println(vsd)

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

    println(id)

    readlines(
        pipeline(
            `docker exec --interactive $id bash -c "/opt/hap.py/bin/hap.py $(joinpath(vvt, basename(vt))) $(joinpath(vvqn, basename(vqn))) -f $(joinpath(vbd, basename(bd))) -r $(joinpath(vre, basename(red))) -o $(joinpath(vouh, "hap.py")) --engine-vcfeval-path $vrt --engine-vcfeval-template $vsd"`,
        ),
    )


    Fastq.support.remove_docker_container(id)

end
