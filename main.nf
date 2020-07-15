nextflow.preview.dsl = 2

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}
params.motif = params.genome ? params.genomes[ params.genome ].motif ?: false : false
params.motif_annot = params.genome ? params.genomes[ params.genome ].motif_annot ?: false : false
params.tss = params.genome ? params.genomes[ params.genome ].tss ?: false : false
params.size = params.genome ? params.genomes[ params.genome ].size ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false

process deg {
    publishDir "$params.outdir/edger", saveAs: { 
        filename -> filename.endsWith(".pdf") 
        ? "plots/$filename" 
        : "lists/$filename"
    }
    input:
    file counts

    output:
    file "*.{pdf,csv}"
    script:
    """
    edger.R $counts
    """
}

process get_deg_promoter {
    input:
    file deg
    output:
    file "*.bed"
    script:
    """
    cut -d, -f1 $deg \
        | tail -n +2 \
        | awk -F, '{print \$1"\t"\$1"_1"}'  \
        | sort -k2,2 \
        | join -1 2 -2 4 - <(sort -k4,4 $params.tss) -t\$'\t' \
        | awk '{OFS="\t";print \$3,\$4-$params.tss_len,\$5+$params.tss_len,\$2,\$6,\$7}' \
        > ${deg.baseName}.bed
    """
}

process get_all_promoter {
    output:
    file "*.bed"
    script:
    """
    cat $params.tss \
        | awk -F"[\t_]" '{OFS="\t";print \$1,\$2-$params.tss_len,\$3+$params.tss_len,\$4,\$6,\$7}' \
        | shuf -n 2000 \
        > ${params.genome}.promoter.bed
    """
}

include fimo from './nf_modules/meme/fimo.nf' params (
    motif: "$params.motif",
    motif_annot: "$params.motif_annot",
    outdir: "$params.outdir"
)

include ame from './nf_modules/meme/ame.nf' params (
    motif: "$params.motif",
    motif_annot: "$params.motif_annot",
    outdir: "$params.outdir"
)

process find_coregulators {
  publishDir "${params.outdir}/coregulators"
  input:
  tuple val(name), file(fimo), val(target)
  output:
  path "*.csv"
  script:
  """
  find_coregulators.py -f $fimo -t $target -a $params.motif_annot -i 10 -o 40 -O ${name}.coregulators_of_${target}.csv
  """
}

include {get_fasta; get_fasta as get_fasta_1} \
    from './nf_modules/bedtools/get_fasta.nf' \
    params (fasta: "$params.fasta")

workflow {
    counts = Channel
        .fromPath( params.counts )
        .ifEmpty { exit 1, "Counts file not found: ${params.counts}" }

    promoter_fa = get_all_promoter | get_fasta

    counts \
        | deg \
        | flatten \
        | filter ( ~/^.*top100_deg.csv/ ) \
        | get_deg_promoter \
        | get_fasta_1 \
        | fimo
    
    get_fasta_1.out \
        | combine (promoter_fa) \
        | ame

    // invariant: co-regulators 
    interested_tfs = Channel.fromPath('./interested_tfs.txt') \
        | splitCsv \
        | map { x -> x[0] } \

    fimo.out \
        | filter ( ~/^.*top.*/ ) \
        | combine (interested_tfs) \
        | find_coregulators

}
