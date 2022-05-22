genome=/home/zekun/lt/moso_genome/Bamboo.HIC.genome.fasta
transcriptome=/home/zekun/lt/moso_genome/bamboo_transcripts.fa
gff=/home/zekun/lt/moso_genome/Bamboo.Hic.gff
gtf=/home/zekun/lt/moso_genome/Bamboo.Hic_gffread_convert.gtf
fastq=/home/zekun/lt/phe_shoot_DRS/fast5_data/
output=/home/zekun/lt/phe_shoot_DRS/nanopore_output/
tools=/home/zekun/lt/tools/
script=/home/zekun/lt/phe_shoot_DRS/phe_shoot_nanopore_script/python/

if [ ! -d "$output/01_fasta" ]; then
 	echo
	echo
	echo "[`date`] Fastq to fasta and fasta(dna) start"
	echo '-----------------------------------------------'
	mkdir $output/01_fasta/
	for i in phe_shoot_2B_rep1 phe_shoot_2M_rep1 phe_shoot_4M_rep1 phe_shoot_4B_rep1
	do
		sed -n '1~4s/^@/>/p;2~4p' $fastq/${i}/fastq/${i}_all.fastq > $output/01_fasta/${i}_all.fa
		python $tools/script/U2T.py $output/01_fasta/${i}_all.fa $output/01_fasta/${i}_all_T.fa
	done
	echo
	echo '-----------------------------------------------'
	echo "[`date`] Fastq to fasta end"
fi 

if [ ! -d "$output/02_correct" ]; then
	echo
	echo
	echo "[`date`] Lordec correct start"
	echo '-----------------------------------------------'
	mkdir $output/02_correct/
	mkdir $output/02_correct/graph_creation_temporary_files
	#Building and saving the de Bruijn Graph of short reads.
	#$tools/lordec-bin_0.9_linux64/lordec-build-SR-graph -T 20 -O $output/02_correct/graph_creation_temporary_files -2 /home/zekun/lt/phe_shoot_illumina/01_rawdata/meta-file.txt -k 19 -s 3 -g $output/02_correct/lordec_SR-graph_k19_s3.h5
	#Correcting the Nanopore reads.
	for i in phe_shoot_2B_rep1 phe_shoot_2M_rep1 phe_shoot_4M_rep1 phe_shoot_4B_rep1
	do 
		$tools/lordec-bin_0.9_linux64/lordec-correct -T 20 -k 19 -s 3 -2 /home/zekun/lt/phe_shoot_illumina/01_rawdata/meta-file.txt -S $output/02_correct/${i}_statistics.log -i $output/01_fasta/${i}_all_T.fa -o $output/02_correct/${i}_lordec.fa -O $output/02_correct/graph_creation_temporary_files > $output/02_correct/${i}_correct.log 2>&1
		rm /home/zekun/lt/phe_shoot_illumina/01_rawdata/meta-file.txt_k19_s3.h5
	done
	echo
	echo '-----------------------------------------------'
	echo "[`date`] Lordec correct end"
fi

if [ ! -d "$output/03_minimap2" ]; then
	echo
	echo
	echo "[`date`] Minimap2 start"
	echo '-----------------------------------------------'
	mkdir $output/03_minimap2/
	for i in phe_shoot_2B_rep1 phe_shoot_2M_rep1 phe_shoot_4M_rep1 phe_shoot_4B_rep1
	do	
		#Alignment
		/home/zekun/lt/tools/minimap2-2.17_x64-linux/minimap2 -ax splice -uf -k14 -t 20 $genome $output/02_correct/${i}_lordec.fa > $output/03_minimap2/${i}_genome_correct.sam
		/home/zekun/lt/tools/minimap2-2.17_x64-linux/minimap2 -ax splice -uf -k14 -t 20 $genome $output/01_fasta/${i}_all_T.fa > $output/03_minimap2/${i}_genome_raw.sam
		/home/zekun/lt/tools/minimap2-2.17_x64-linux/minimap2 -ax map-ont -uf -t 20 $transcriptome $output/02_correct/${i}_lordec.fa > $output/03_minimap2/${i}_transcriptome_correct.sam
		/home/zekun/lt/tools/minimap2-2.17_x64-linux/minimap2 -ax map-ont -uf -t 20 $transcriptome $output/01_fasta/${i}_all_T.fa > $output/03_minimap2/${i}_transcriptome_raw.sam
		
		#SAM to BAM
		samtools view -Sb -o $output/03_minimap2/${i}_genome_correct.bam $output/03_minimap2/${i}_genome_correct.sam
		samtools view -Sb -o $output/03_minimap2/${i}_genome_raw.bam $output/03_minimap2/${i}_genome_raw.sam
		samtools view -Sb -o $output/03_minimap2/${i}_transcriptome_correct.bam $output/03_minimap2/${i}_transcriptome_correct.sam
		samtools view -Sb -o $output/03_minimap2/${i}_transcriptome_raw.bam $output/03_minimap2/${i}_transcriptome_raw.sam

		#Sort BAM
		samtools sort $output/03_minimap2/${i}_genome_correct.bam $output/03_minimap2/${i}_genome_correct.sorted.bam
		samtools sort $output/03_minimap2/${i}_genome_raw.bam $output/03_minimap2/${i}_genome_raw.sorted.bam
		samtools sort $output/03_minimap2/${i}_transcriptome_correct.bam $output/03_minimap2/${i}_transcriptome_correct.sorted.bam
		samtools sort $output/03_minimap2/${i}_transcriptome_raw.bam $output/03_minimap2/${i}_transcriptome_raw.sorted.bam

		#Index
		samtools index $output/03_minimap2/${i}_genome_correct.sorted.bam
		samtools index $output/03_minimap2/${i}_genome_raw.sorted.bam
		samtools index $output/03_minimap2/${i}_transcriptome_correct.sorted.bam
		samtools index $output/03_minimap2/${i}_transcriptome_raw.sorted.bam
	done
	echo
	echo '-----------------------------------------------'
	echo "[`date`] Minimap2 end"
fi

if [ ! -d "$output/04_stringtie" ]; then
	echo
	echo
	echo "[`date`] Stringtie start"
	echo '-----------------------------------------------'
	mkdir $output/04_stringtie/
	for i in phe_shoot_2B_rep1 phe_shoot_2M_rep1 phe_shoot_4M_rep1 phe_shoot_4B_rep1
	do	
		#-F 2308. Discard unmapped/supplementary_alignment/not_primary_alignment reads.
	  	samtools view -Sb -F 2308 -o $output/03_minimap2/${i}_genome_correct_2308.bam $output/03_minimap2/${i}_genome_correct.sam
		samtools sort $output/03_minimap2/${i}_genome_correct_2308.bam $output/03_minimap2/${i}_genome_correct_2308.sorted.bam
		samtools index $output/03_minimap2/${i}_genome_correct_2308.sorted.bam
		rm $output/03_minimap2/${i}_genome_correct_2308.bam

		$tools/stringtie-2.1.4.Linux_x86_64/stringtie -L -G $gff -o $output/04_stringtie/${i}_stringtie.gtf $output/03_minimap2/${i}_genome_correct_2308.sorted.bam
	done
		
	#Merge gtf file
	ls $output/04_stringtie/*.gtf > $output/04_stringtie/stringtie_gff_list.txt
	$tools/stringtie-2.1.4.Linux_x86_64/stringtie --merge -o $output/04_stringtie/phe_nanopore_merge.gtf -G -p 10 $output/04_stringtie/stringtie_gff_list.txt	
	echo
	echo '-----------------------------------------------'
	echo "[`date`] Stringtie end"
fi

if [ ! -d "$output/05_AS" ]; then
	echo
	echo
	echo "[`date`] AS(SUPPA) detection start"
	echo '-----------------------------------------------'
	mkdir $output/05_AS/
	source /home/zekun/anaconda3/bin/activate suppa
	for i in phe_shoot_2B_rep1 phe_shoot_2M_rep1 phe_shoot_4M_rep1 phe_shoot_4B_rep1 
	do
		python /home/zekun/anaconda3/envs/suppa/bin/suppa.py generateEvents -i $output/04_stringtie/${i}_stringtie.gtf -o $output/05_AS/ -f ioe -e SE SS MX RI FL
	done
	python /home/zekun/anaconda3/envs/suppa/bin/suppa.py generateEvents -i $output/04_stringtie/phe_nanopore_merge.gtf -o $output/05_AS/ -f ioe -e SE SS MX RI FL
	conda deactivate
	echo
	echo '-----------------------------------------------'
	echo "[`date`] AS(SUPPA) detection end"
fi

if [ ! -d "$output/06_featureCounts" ]; then
	echo
	echo
	echo "[`date`] featureCounts start"
	echo '-----------------------------------------------'
	mkdir $output/06_featureCounts/
	for i in phe_shoot_2B_rep1 phe_shoot_2M_rep1 phe_shoot_4M_rep1 phe_shoot_4B_rep1
        do
		$tools/subread-2.0.1-Linux-x86_64/bin/featureCounts -a $gtf -L -R CORE -t exon -g gene_id -o $output/06_featureCounts/${i}_counts.txt $output/03_minimap2/${i}_genome_correct_2308.sorted.bam > $output/06_featureCounts/${i}_counts.log 2>&
	done
	echo
	echo '-----------------------------------------------'
	echo "[`date`] featureCounts end"
fi

if [ ! -d "$output/07_readsSTAT" ]; then
	echo
	echo
	echo "[`date`] 07_readsSTAT start"
	echo '-----------------------------------------------'
	mkdir $output/07_readsSTAT/
	for i in phe_shoot_2B_rep1 phe_shoot_2M_rep1 phe_shoot_4M_rep1 phe_shoot_4B_rep1
	do
		python $script/reads_coverage.py $output/03_minimap2/${i}_transcriptome_correct.sorted.bam $output/07_readsSTAT/{i}
	done
	#Statistics of read length
	python $script/reads_statistics.py $output/03_minimap2/phe_shoot_2B_rep1_transcriptome_correct.sorted.bam,$output/03_minimap2/phe_shoot_2M_rep1_transcriptome_correct.sorted.bam,$output/03_minimap2/phe_shoot_4B_rep1_transcriptome_correct.sorted.bam,$output/03_minimap2/phe_shoot_4M_rep1_transcriptome_correct.sorted.bam phe_shoot_2B_rep1,phe_shoot_2M_rep1,phe_shoot_4B_rep1,phe_shoot_4M_rep1 > $output/07_readsSTAT/nanopore_correct_reads.stat
	python $script/reads_statistics.py $output/03_minimap2/phe_shoot_2B_rep1_transcriptome_raw.sorted.bam,$output/03_minimap2/phe_shoot_2M_rep1_transcriptome_raw.sorted.bam,$output/03_minimap2/phe_shoot_4B_rep1_transcriptome_raw.sorted.bam,$output/03_minimap2/phe_shoot_4M_rep1_transcriptome_raw.sorted.bam phe_shoot_2B_rep1,phe_shoot_2M_rep1,phe_shoot_4M_rep1,phe_shoot_4B_rep1 > $output/07_readsSTAT/nanopore_raw_reads.stat
	echo
	echo '-----------------------------------------------'
	echo "[`date`] 07_readsSTAT end"
fi

if [ ! -d "$output/08_bedtools_coverage" ]; then
	echo
	echo
	echo "[`date`] 08_bedtools_coverage start"
	echo '-----------------------------------------------'
	mkdir $output/08_bedtools_coverage/

	for i in phe_shoot_2B_rep1 phe_shoot_2M_rep1 phe_shoot_4M_rep1 phe_shoot_4B_rep1
	do
		samtools view -Sb -F 2308 -o $output/03_minimap2/${i}_transcriptome_correct_2308.bam $output/03_minimap2/${i}_transcriptome_correct.sam
		samtools sort $output/03_minimap2/${i}_transcriptome_correct_2308.bam $output/03_minimap2/${i}_transcriptome_correct_2308.sorted.bam
		samtools index $output/03_minimap2/${i}_transcriptome_correct_2308.sorted.bam
		rm $output/03_minimap2/${i}_transcriptome_correct_2308.bam
	done
	# Merge bam file.
	samtools merge $output/03_minimap2/phe_shoot_all_rep1_transcriptome_correct_2308.sorted.bam $output/03_minimap2/phe_shoot_2B_rep1_transcriptome_correct_2308.sorted.bam,$output/03_minimap2/phe_shoot_2M_rep1_transcriptome_correct_2308.sorted.bam,$output/03_minimap2/phe_shoot_4B_rep1_transcriptome_correct_2308.sorted.bam,$output/03_minimap2/phe_shoot_4M_rep1_transcriptome_correct_2308.sorted.bam

	source /home/zekun/anaconda3/bin/activate nanom6a
	for i in phe_shoot_2B_rep1 phe_shoot_2M_rep1 phe_shoot_4M_rep1 phe_shoot_4B_rep1 phe_shoot_all_rep1
	do
		bedtools coverage -d -split -a /home/zekun/lt/moso_genome/Bamboo.Hic.transcritome.bed -b $output/03_minimap2/${i}_transcriptome_correct_2308.sorted.bam > $output/08_bedtools_coverage/${i}_bedtools_coverage.bed
	done
	conda deactivate
	echo
	echo '-----------------------------------------------'
	echo "[`date`] 08_bedtools_coverage end"
fi
















