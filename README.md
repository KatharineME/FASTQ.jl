# Fastq.jl

Turns raw sequencing reads into intpretable files like: VCF for variants or a gene by sample matrix for gene expression analysis.


For each of the use cases below, there is a corresponding Jupyter Notebook workflow:

Germline DNA > VCF (Genome alignment)

Somatic DNA > VCF (Genome alignment)

cDNA > VCF (Genome alignment)

cDNA > Gene by sample (Psuedoalignment)


## Use

```jl
using Fastq.jl
```

### `function1`

```jl
align()
```

Minimap2 flags

`-x ` 
A recommended meta flag that specifies alginment chain bandwidth, elongation, discard, how mnay secondary alignment should be output, and more. 
`-a` 
Generates CIGAR and output alignments in SAM format.
`--sr` 
Enable short read mode.
`--splice` Enable splice alignment mode.
`-uf`
Use transcript strand to find canonical splice sites. 
`-t NUM`
Number of threads
`-K NUM`
Number of bases loaded into memory to process in a mini-batch [500M]. A large NUM pased here helps with load balancing in the multi-threading mode.
`-R`
SAM read group line in a format like @RG\\tID:foo\\tSM:bar [].

samtools fixmate

## Contribution

To report a bug, request a feature, or leave a comment (about anything related to this repository), just [submit an issue](https://github.com/KatharineME/Fastq.jl.jl/issues/new/choose).

---

Made by https://github.com/KwatMDPhD/PkgRepository.jl
