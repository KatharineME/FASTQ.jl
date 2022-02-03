function call_somatic_variant(
    mo::String,
    ge::String,
    so::String,
    fa::String,
    chs::String,
    chn::String,
    pao::String,
    n_jo::Int,
    me::Int,
    to::String,
    pas::String,
)::Nothing

    if !(isfile("$fa.fai") && ispath("$fa.gzi"))

        run(`samtools faidx $fa`)

    end
    
    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end
   
    if check_directory(pao, "call somatic variant")

        return nothing

    end


    # Run docker container
    
    id, voo, vof, voc, voge, vot = run_docker_container(to, fa, chs, ge, pao)


    println("docker is running")


    # Set config parameters

    co = "--referenceFasta $fa --callRegions $chs --normalBam $ge --tumorBam $so"


    # Set run parameters

    ruString = "--mode local --jobs $n_jo --memGb $me --quiet"

    pav = joinpath("results", "variants")



    # Configure and run manta

    pam = joinpath(pao, "manta")

    pamr = joinpath(pam, "runWorkflow.py")

    sc = "manta-1.6.0.centos6_x86_64/bin/configManta.py" 

    se =  readlines(pipeline(`docker exec --interactive $id bash -c "./home/$vot/$(sc) $co --outputContig --runDir /home/$pam && ./home/$pamr $ru"`))

    println("$(join(re, " "))\n")

    println("configured docker")



    # Configure and run strelka

    pas = joinpath(pao, "strelka")
    
    st = "configureStrelkaSomaticWorkflow.py $co --indelCandidates $(joinpath(pam, pav, "candidateSmallIndels.vcf.gz")) --runDir $pas"

    pasr = joinpath(pas, "runWorkflow.py")

    run(`bash -c "source activate py2 && $st && $pasr $ru"`)



    # bcftools

    sa = joinpath(pao, "sample.txt")

    # TODO: get sample names (maybe from .bam) and use them instead of "Germ" and "Soma"

    open(io -> write(io, "Germ\nSoma"), sa; write = true)

    pain = joinpath(pas, pav, "somatic.indels.vcf.gz")

    run(
        pipeline(
            `bcftools reheader --threads $n_jo --samples $sa $pain`,
            "$pain.tmp",
        ),
    )

    mv("$pain.tmp", pain; force = true)

    run(`tabix --force $pain`)

    pasv = joinpath(pas, pav, "somatic.snvs.vcf.gz")

    run(
        pipeline(
            `bcftools reheader --threads $n_jo --samples $sa $pasv`,
            "pasv.tmp",
        ),
    )

    mv("$pasv.tmp", pasv; force = true)

    run(`tabix --force $pasv`)

    vc_ = [joinpath(pam, pav, "somaticSV.vcf.gz"), pain, pasv]



    paco = joinpath(pao, "concat.vcf.gz")

    run(
        pipeline(
            `bcftools concat --threads $n_jo --allow-overlaps $vc_`,
            `bcftools annotate --threads $n_jo --rename-chrs $chn`,
            `bgzip --threads $n_jo --stdout`,
            paco,
        ),
    )

    run(`tabix $paco`)

    

    # snpeff

    sn = joinpath(pao, "snpeff")

    mkpath(sn)

    snvc = joinpath(sn, "snpeff.vcf.gz")

    run(
        pipeline(
            `java -Xmx$(me)g -jar $pas GRCh38.99 -noLog -verbose -csvStats $(joinpath(sn, "stats.csv")) -htmlStats $(joinpath(sn, "stats.html")) $paco`,
            `bgzip --threads $n_jo --stdout`,
            snvc,
        ),
    )

    run(`tabix $snvc`)

    ps = joinpath(pao, "pass.vcf.gz")

    run(
        pipeline(
            `bcftools view --threads $n_jo --include 'FILTER=="PASS"' $snvc`,
            `bgzip --threads $n_jo --stdout`,
            ps,
        ),
    )

    return run(`tabix $ps`)

end
