#!/usr/env python3

import sys
import subprocess

DELETE_INTERMEDIATES = 'yes'

INPUT_FOLDER = '/Volumes/Cocoon/Killifish_ATACseq_Pipeline/Input'
INTERMEDIATE_FOLDER = '/Volumes/Cocoon/Killifish_ATACseq_Pipeline/Intermediate'
OUTPUT_FOLDER = '/Volumes/Cocoon/Killifish_ATACseq_Pipeline/Output'
FIGURE_FOLDER = '/Volumes/Cocoon/Killifish_ATACseq_Pipeline/Figure'
TSS_FOLDER = '/Volumes/Cocoon/Killifish_ATACseq_Pipeline/Input/TSS'
PICARD_JAR = '/Volumes/Cocoon/Killifish_ATACseq_Pipeline/picard.jar'

GENOME = 'Nfur.fa'
FASTQ1_BASE = 'Nfur_1M_1_ATAC'
FASTQ2_BASE = 'Nfur_1M_2_ATAC'
FASTQ3_BASE = 'Nfur_preDia_1_ATAC'
FASTQ4_BASE = 'Nfur_preDia_2_ATAC'
BASES = [FASTQ1_BASE, FASTQ2_BASE, FASTQ3_BASE, FASTQ4_BASE]

'''
#Create genome index for read-mapping to the genome
subprocess.call(['bowtie2-build', '-f', INPUT_FOLDER + '/' + GENOME,  INPUT_FOLDER + '/NFUR'])


#Trim reads from paired fastq files and generate fastqc reports of read quality
for item in BASES:
	FQCA = '\"--outdir ' + FIGURE_FOLDER + '/\"'
	subprocess.call(['trim_galore', '--fastqc', '--fastqc_args', FQCA, '--gzip', '--output_dir', INPUT_FOLDER + '/',  '--paired', INPUT_FOLDER + '/' + item + '_R1.fastq.gz', INPUT_FOLDER + '/' + item + '_R2.fastq.gz'])


#Map reads from each fastq file pair to the genome and log the alignment rates
with open(OUTPUT_FOLDER + '/align_log.txt', 'wt') as logger:
	for item in BASES:
		logger.write(item + '\n')
		subprocess.call(['bowtie2', '--very-sensitive', '-x', INPUT_FOLDER + '/NFUR', '-1', INPUT_FOLDER + '/' + item + '_R1_val_1.fq.gz', '-2', INPUT_FOLDER + '/' + item + '_R2_val_2.fq.gz', '-S', INTERMEDIATE_FOLDER + '/' + item + '.sam'], stdout = logger)

#Convert aligned read SAM files to BAM format and sort the reads by chromosome
for item in BASES:
	with open(INTERMEDIATE_FOLDER + '/' + item + '.bam', 'wt') as bammer:
		subprocess.call(['samtools', 'view', '-S', '-b', INTERMEDIATE_FOLDER + '/' + item + '.sam'], stdout=bammer) 
		subprocess.call(['samtools', 'sort', '-T', INTERMEDIATE_FOLDER + '/' + item + 'temp.bam', '-o', INTERMEDIATE_FOLDER + '/' + item + '_sorted.bam',  INTERMEDIATE_FOLDER + '/' + item + '.bam'])
		subprocess.call(['samtools', 'index', INTERMEDIATE_FOLDER + '/' + item + '_sorted.bam'])

#Collect metrics on the presence of mitochondrial mapped reads and remove them from the BAM files
with open(OUTPUT_FOLDER + '/mito_read_percent.log', 'wt') as logger3:
	logger3.write('Sample\tMitochondrial_Read_Percentage\n')
	
	for item in BASES:
	
		with open(INTERMEDIATE_FOLDER + '/' + item + '_mitotemp.sam', 'wt') as MT, open(INTERMEDIATE_FOLDER + '/' + item + '_alltemp.sam', 'wt') as AT:
			subprocess.call(['samtools', 'view', '-c', INTERMEDIATE_FOLDER + '/' + item + '_sorted.bam', 'NC_011814.1'], stdout = MT)
			subprocess.call(['samtools', 'view', '-c', INTERMEDIATE_FOLDER + '/' + item + '_sorted.bam'], stdout = AT)
			
		with open(INTERMEDIATE_FOLDER + '/' + item + '_mitotemp.sam', 'r') as MT, open(INTERMEDIATE_FOLDER + '/' + item + '_alltemp.sam', 'r') as AT:
			MT_Count = len(MT.readlines())
			AT_Count = len(AT.readlines())
			MT_Percent = 100*(MT_Count/AT_Count)
			logger3.write(item + '\t' + str(MT_Percent) + '\n')
		
		REMOVER = 'samtools view -h ' + INTERMEDIATE_FOLDER + '/' + item + '_sorted.bam' + ' | grep -v NC_011814.1 | samtools sort -O bam -o ' + INTERMEDIATE_FOLDER + '/' + item + '_sorted_mitoremoved.bam'
		subprocess.Popen(['sh', '-c', REMOVER], stdout=subprocess.PIPE) 		
'''

#Use picard tools to mark all duplicate reads and collect metrics on them
for item in BASES:
	subprocess.call(['java', '-jar', PICARD_JAR, 'MarkDuplicates', 'QUIET=true', 'INPUT=' + INTERMEDIATE_FOLDER + '/' + item + '_sorted_mitoremoved.bam', 'OUTPUT=' + INTERMEDIATE_FOLDER + '/' + item + '_sorted_mitoremoved_dupmarked.bam', 'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000', 'METRICS_FILE=' + OUTPUT_FOLDER + '/' + item + 'dupmarked.metrics', 'REMOVE_DUPLICATES=false', 'CREATE_INDEX=true', 'VALIDATION_STRINGENCY=LENIENT', 'TMP_DIR=' + INTERMEDIATE_FOLDER + '/tmp/'])

'''
#Filter duplicated reads and reads with suboptimal quality, resort the reads, and index the file for visualization
for item in BASES:

	with open(INTERMEDIATE_FOLDER + '/' + item + 'mitoremoved_dupremoved.bam', 'wt') as dupper: 
		subprocess.call(['samtools', 'view', '-h', '-b', '-F', '1024', INTERMEDIATE_FOLDER + '/' + item + '_sorted_mitoremoved_dupmarked.bam'], stdout = dupper)
		
	with open(INTERMEDIATE_FOLDER + '/' + item + '_mitoremoved_dupremoved_mapq20.bam', 'wt') as mapper:
		subprocess.call(['samtools', 'view', '-h', '-q', '20', INTERMEDIATE_FOLDER + '/' + item + '_mitoremoved_dupremoved.bam'], stdout = mapper)

	with open(INTERMEDIATE_FOLDER + '/' + item + '_filtered.bam', 'wt') as filter:
		subprocess.call(['samtools', 'view', '-h', '-b', '-F', '1804', '-f', '2', INTERMEDIATE_FOLDER + '/' + item + '_mitoremoved_dupremoved_mapq20.bam'], stdout = filter)

	subprocess.call(['samtools', 'sort', '-T', INTERMEDIATE_FOLDER + '/' + item + '_filtered_temp.bam', '-o', INTERMEDIATE_FOLDER + '/' + item + '_filtered_sorted.bam'])


#Shift read alignments to account for Tn5 Binding bias and resort the BAM files
for item in BASES:
	subprocess.call(['alignmentSieve', '--numberOfProcessors', '4', '--ATACshift', '--bam', INTERMEDIATE_FOLDER + '/' + item + '_filtered_sorted.bam', '-o', INTERMEDIATE_FOLDER + '/' + item + '_filtered_shifted.bam'])
	subprocess.call(['samtools', 'sort', '-O', 'bam', '-o', OUTPUT_FOLDER + '/' + item + '_FINAL.bam', INTERMEDIATE_FOLDER + '/' + item + '_filtered_shifted.bam'])
	subprocess.call(['samtools', 'index', OUTPUT_FOLDER + '/' + item + '_FINAL.bam'])


#Call peaks from finalized BAM files
with open(OUTPUT_FOLDER + '/macs2.log', 'wt') as logger2:
	for item in BASES:
		subprocess.call(['macs2', 'callpeak', '-f', 'BAMPE', '-g', '8.568x10^8', '--keep-dup', 'all', '--cutoff-analysis', '-n', 'killifish', '-B', '-t', OUTPUT_FOLDER + '/' + item + '_FINAL.bam', '--outdir', OUTPUT_FOLDER + '/'], stdout = logger2)


#Use the called peaks and feature data to evaluate TSS enrichment
for item in BASES:
	subprocess.call(['computeMatrix', 'reference-point', '--referencePoint', 'TSS', '-b', '2000', '-a', '2000', '-R', INPUT_FOLDER + '/TSS/', '-S', OUTPUT_FOLDER + '/' + item + '_FINAL.bw', '--missingDataAsZero', '--skipZeros', '-o', OUTPUT_FOLDER + '/' + item + '_FINAL_matrix.gz'])
	subprocess.call(['plotHeatmap', '-m', OUTPUT_FOLDER + '/' + item + '_FINAL_matrix.gz', '-out', FIGURE_FOLDER + '/' + item + '_FINAL_matrix.pdf', '--colorList=\'white,blue\'', '--plotFileFormat', 'pdf'])


#Run first R script to get concensus and differential peaks
subprocess.call(['rscript', './Killifish_ATACseq_Pipeline.R', 'RUN1', INPUT_FOLDER, OUTPUT_FOLDER, FIGURE_FOLDER])


#Generate bed files from differential peak sets
######


#Get motifs
subprocess.call(['findMotifsGenome.pl', 'Nfurzeri', OUTPUT_FOLDER + '/DIA_UP_DE.bed', OUTPUT_FOLDER + '/', '-mset', 'vertebrates'])
subprocess.call(['findMotifsGenome.pl', 'Nfurzeri', OUTPUT_FOLDER + '/DEV_UP_DE.bed', OUTPUT_FOLDER + '/', '-mset', 'vertebrates'])


#Organize Motifs
in = knownreuslts.txt
out = knonresults.csv

#Make Motif Heatmap
subprocess.call(['rscript', './Killifish_ATACseq_Pipeline.R', 'RUN2', INPUT_FOLDER, OUTPUT_FOLDER, FIGURE_FOLDER])


#Remove intermediate files
if DELETE_INTERMEDIATES == 'yes':
	subprocess.call(['rm', '-rf', INTERMEDIATE_FOLDER])

'''	