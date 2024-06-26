module HLAType

using ..FASTQ

function hlatype(ba)

    # copy bam to HLA directory


    # run xhla
    # docker run -v `pwd`:`pwd` -w `pwd` humanlongevity/hla --sample_id 1020 --input_bam_path tests/1020.bam --output_path test

    # copy results over to callgermlinevariants

    # delete copied bam file

end

end
