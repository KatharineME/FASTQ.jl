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

    #Fastq.support.error_if_directory(pa)


    # Make vcfeval SDF

    red = strip(re, ".gz")

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

    ## run vcfeval

    ouv = joinpath(pa, "vcfeval")

    run(`$rt vcfeval 
        --baseline=$vt 
        --bed-regions=$bd 
        --calls=$vqn 
        --template=$sd 
        --output=$ouv
        --threads=$n_jo`)


    ## run hap.py in hap.py docker container

    ouh = joinpath(pa, "happy/")

    mkdir(ouh)

    ho = "/home/"

    pvt = joinpath(splitpath(vt)[end - 1], splitpath(vt)[end])

    vvt = joinpath(ho, pvt)

    pvqn = joinpath(splitpath(vqn)[end -1], splitpath(vqn)[end])

    vvqn = joinpath(ho, pvqn)

    vbd = joinpath(ho, basename(bd))

    pre = joinpath(splitpath(red)[end - 1], splitpath(red)[end])

    vre = joinpath(ho, pre)

    vouh = joinpath(ho, splitpath(ouh)[end])

    vrt = joinpath(ho, basename(rt))

    vsd = joinpath(ho, basename(sd))

    println(vvt)

    println(vvqn)

    println(vbd)

    println(vre)

    println(vouh)

    println(vrt)

    println(vsd)

    id = readlines(pipeline(`docker run 
            --interactive 
            --detach 
            --tty 
            --user root
            -v $pvt:$vvt 
            -v $pvqn:$vvqn 
            -v $bd:$vbd 
            -v $pre:$vre 
            -v $ouh:$vouh 
            -v $rt:$vrt 
            -v $sd:$vsd 
            pkrusche/hap.py
            bash`))

    println(id)

    readlines(
        pipeline(
            `docker exec --interactive $id bash -c "/opt/hap.py/bin/hap.py $vvt $vvqn -f $vbd -r $vre -o $(joinpath(vouh, "hap.py")) --engine-vcfeval-path $vrt --engine-vcfeval-template $vsd"`,
        ),
    )

end
