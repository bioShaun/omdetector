#!/usr/bin/env nextflow


def helpMessage() {
    log.info """

    A pipeline to identify unmapped reads source species

    Usage:

    ==================================================================

    References
      --uniprot_db                  Path to uniprot database
      --unmapped_fa                 Path to unmapped fasta
    
    Mandatory arguments:
      --chunk_size                  sequence per file for blast
      --outdir                      Path to analysis results 
    ===================================================================
    """.stripIndent()
}
/*
 * SET UP CONFIGURATION VARIABLES
 */

 
// workflow internal path&files
script_dir = file("$baseDir/script/")

 // Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

def check_file_exist = {file_path, file_type ->
    if (file_path) {
        file_path = file(file_path)
        if( !file_path.exists() ) exit 1, "${file_type} file not found: ${file_path}"
        return file_path
    } else {
        exit 1, "${file_type} is required!"
    }
}

// load parameters
params.uniprot_db = '/public/database/swissprot/uniprot_sprot.fasta'
params.chunk_size = 100
params.unmapped_fa = false
params.outdir = false

db_name = file(params.uniprot_db).name
db_path = file(params.uniprot_db).parent

input_fa = check_file_exist(params.unmapped_fa, 'input fasta')

Channel
    .fromPath(input_fa)
    .splitFasta(by: params.chunk_size)
    .set { split_fa }


// orf predict
process orf_predict {

    input:
    file fa from split_fa

    output:
    file "${fa}.transdecoder_dir/longest_orfs.pep" into orf_fa

    cpus = 4

    script:
    """
    TransDecoder.LongOrfs -t ${fa}
    """
}



// run blast
process run_blast {

    input:
    file orf_fa from orf_fa
    file db_path 
 
    output:
    file "${orf_fa}.blastout" into blast_result

    cpus = 40
 
    """
    blastp -db ${db_path}/${db_name} -query ${orf_fa} \\
        -evalue 1e-5 \\
        -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen evalue bitscore stitle" \\
        -max_target_seqs 1 -num_threads ${task.cpus} > ${orf_fa}.blastout
    """
}  

// merge blastout
outdir = file(params.outdir)
outdir.mkdirs()
blastout = file("${params.outdir}/merged.blastout")

blast_result
    .collectFile(name: blastout)
    .set{ m_blast_result }

// cal stat
process seq_sp_info {

    publishDir "${outdir}", mode: 'copy'

    input:
    file blastout from m_blast_result

    output:
    file 'species.stat.txt'
 
    """
    python ${script_dir}/sp_stat.py --blastout ${blastout} --outfile species.stat.txt
    """
}  
