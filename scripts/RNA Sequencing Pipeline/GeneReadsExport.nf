#!/usr/bin/env nextflow


params.fastqs = "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/Genomics_Innovation_Unit/internship_proper/sra_data/ecoli/SRR139*_{1,2}.fastq"
params.gtf = "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/Genomics_Innovation_Unit/internship_proper/sra_data/ecoli/GCF_003017915.1_ASM301791v1_genomic.gtf"
params.reference = "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/Genomics_Innovation_Unit/internship_proper/sra_data/ecoli/GCF_003017915.1_ASM301791v1_genomic.fna"
params.outdir = "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/Genomics_Innovation_Unit/internship_proper/sra_data/ecoli"
params.genes = "/Users/rugaremaruzani/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/Genomics_Innovation_Unit/internship_proper/sra_data/ecoli/gene_file.txt"

Channel.fromPath(params.reference).set { ref_genome_ch }

reads_pairs_ch = Channel.fromFilePairs(params.fastqs)


process STAR_align {
       input:
       tuple pair_id, path(fastqs) from reads_pairs_ch

       output:
       path '*.sam' into output_ch

       script:
       """
       STAR --runThreadN $task.cpu \
               --genomeDir $params.outdir \
               --readFilesIn ${fastqs[0]} ${fastqs[1]}
       
       mv 'Aligned.out.sam' '${pair_id}'.sam
       """
}


process sam_bam {
       input:
       path sam_file from output_ch

	output:
	path "*.bam" into bam_ch

       script:
       """
        samtools view -bS -o "${sam_file.baseName}.bam" ${sam_file}
       """
}

process gene_gff {
       input:
       path gtf from params.gtf
	path file from params.genes

        output:
        path ("*.gff") into (gene_gff_ch1, gene_gff_ch2)

       shell:
       """
	while read -r line; do grep -w \${line}  $params.gtf > "\${line}.gff"; done < gene_file.txt
       """
}

process gene_bed {
       input:
       each gffs from gene_gff_ch1
          
        output:
        path ("*.bed") into gene_bed_ch

       shell:
       """ 
	agat_convert_sp_gff2bed.pl --gff ${gffs} -o  ${gffs.baseName}.bed
       """ 
}



process gene_bam {

       input:
       each beds from gene_bed_ch
       val bam_file from bam_ch
       each gff from gene_gff_ch2


        output:
        path ("*.bam") into (bam_ch1, bam_ch2)

       shell:
       """
       samtools view -b -L ${beds} -o ${beds.baseName}.${bam_file.baseName}.bam $bam_file
	"""
}

process export_reads {
	publishDir "${params.outdir}/gene_fastqs", mode: 'copy'

       input:
	each bam from bam_ch1

        output:
	tuple val(bam), path("*.fastq")
	
       shell:
       """
        gatk SamToFastq -INTER true -I ${bam} -F ${bam.baseName}.fastq --VALIDATION_STRINGENCY SILENT
       """
}



process sanity_check {
       publishDir "${params.outdir}/sanity_check", mode: 'copy'

       input:
	each bam from bam_ch2

        output:
        path ("*.txt") into sanity_ch

       shell:
       """
        featureCounts --extraAttributes gene -p -T 5 -t gene -a $params.gtf -o ${bam.baseName}.counts.txt ${bam}
	grep -wv "0" ${bam.baseName}.counts.txt > ${bam.baseName}.count.txt	
	rm ${bam.baseName}.counts.txt
       """
}



