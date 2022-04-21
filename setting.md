`number_of_jobs` _Int_<br>
Number of jobs

`memory` _Int_<br>
Gigabytes of memory

`organism` _String_<br>
"human" or "mouse" used by `apply_cdna_to_transcriptome` when creating gene_x_sample matrix

`molecule` _String_<br>
"dna" or "cdna"

`exome` _Bool_<br>
"true" or "false"

`fragment_length` _Int_<br>
The most common cDNA fragment length in your data, used by `apply_cdna_to_genome`

`fragment_length_standard_deviation` _Int_<br>
For cDNA alignment, the estimated standard deviation of fragment lengths

`sample` _String_<br>
Sample name used for file naming and name of sample column in VCF

`annotate_with_rsid` _Bool_<br>
"true" or "false"

`tool_directory` _String_<br>
Path to directory storing strelka and manta

`output_directory` _String_<br>
Path to where you want output files

`read_name_scheme` _String_<br>
Naming scheme your reads use, "R1" or "read1" for example

`read1` _String_<br>
Path to DNA forward read file used by `apply_germline_dna_to_genome` and `apply_somatic_dna_to_genome`

`read2` _String_<br> 
Path to DNA reverse read file used by `apply_germline_dna_to_genome` and `apply_somatic_dna_to_genome`

`somatic_read1` _String_<br>
Path to somatic DNA forward read file used by `apply_somatic_dna_to_genome`

`somatic_read2` _String_<br>
Path to somatic DNA reverse read file used by `apply_somatic_dna_to_genome`

`dna_read_directory` _String_<br>
Path to directory with DNA read files used by `concatenate_fastq` to combine fastqs of the same read direction

`cdna_read_directory` _String_<br>
Path to directory with cDNA read files used by `apply_cdna_to_genome` and `apply_cdna_to_transcriptome` which both use the read name and directory structure to keep different cDNA samples separate. 

`reference_genome` _String_<br>
Path to human reference genome fna.gz or fa.gz file used by `apply_germline_dna_to_genome`, `apply_somatic_dna_to_genome`, and `apply_cdna_to_genome`

`reference_transcriptome` _String_<br>
Path to human reference transcriptome fna.gz or fa.gz file used by `apply_cdna_to_transcriptome`

`mouse_transcript_to_mouse_gene` _String_<br>
Path to tsv with mouse transcript name column and mouse gene name column used by `apply_cdna_to_transcriptome` when `organism` is set to "mouse"

`chromosome_position` _String_<br>
Path to tsv with 3 columns: chromosome name, start position of chromosome (0), and end position of chromosome. Used by `apply_germline_dna_to_genome`, `apply_somatic_dna_to_genome`, and `apply_cdna_to_genome`

`chromosome_name` _String_<br>
Path to tsv with column of chromosome names ("chr1", "chr2") and column of their integer names (1, 2)

`variant_database` _String_<br>
Path to variant database VCF file (dbsnp or ensembl variant database for example) used by `apply_germline_dna_to_genome`, `apply_somatic_dna_to_genome`, and `apply_cdna_to_genome` when `annotate_with_rsid` is set to "true"

`snpeff` _String_<br>
Path to snpEff.jar file used by `apply_germline_dna_to_genome`, `apply_somatic_dna_to_genome`, and `apply_cdna_to_genome`
