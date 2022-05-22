input=/home/lfgu/bigdata/hitseq/phe_shoot_nanopore/Illumina_raw_data/
tools=/home/lfgu/bigdata/hitseq/phe_shoot_nanopore/tool/
output=/home/lfgu/bigdata/hitseq/phe_shoot_nanopore/illumina_output/
samtools=/home/lfgu/bigdata/hitseq/phe_shoot_nanopore/tool/miniconda3/envs/samtools/bin/
base_bin=/home/lfgu/bigdata/hitseq/phe_shoot_nanopore/tool/miniconda3/bin/
hisat2_index=/home/lfgu/bigdata/hitseq/phe_shoot_nanopore/input/Bamboo.HIC.genome.HISAT2.index/Bamboo.HIC.genome
bamboo_gtf=/home/lfgu/bigdata/hitseq/phe_shoot_nanopore/input/Bamboo.Hic_gffread_convert.gtf
bamboo_gff=/home/lfgu/bigdata/hitseq/phe_shoot_nanopore/input/Bamboo.Hic.gff
featurecounts=/home/lfgu/bigdata/hitseq/phe_shoot_nanopore/tool/subread-2.0.1-Linux-x86_64/bin/
susupplementary=/home/lfgu/bigdata/hitseq/phe_shoot_nanopore/illumina_output/supplementary/

if [ ! -d "$output/01_fastp/" ]; then
        echo
        echo
        echo "[`date`] fastq quality control"
        echo '-----------------------------------------------'
        mkdir $output/01_fastp/
	mkdir $output/01_fastp/report/
	mkdir $output/01_fastp/log/
    	
	$base_bin/parallel --jobs 8 $tools/fastp --in1 $input/Phe-{1}-Rep1/Phe-{1}-Rep1_{2}.R1.fastq --in2 $input/Phe-{1}-Rep1/Phe-{1}-Rep1_{2}.R2.fastq --out1 $output/01_fastp/Phe-{1}-Rep1_{2}.R1.fastp.fastq --out2 $output/01_fastp/Phe-{1}-Rep1_{2}.R2.fastp.fastq -h $output/01_fastp/report/Phe-{1}-Rep1_{2}_fastp.html --json $output/01_fastp/report/Phe-{1}-Rep1_{2}_fastp.json --thread 5 '>' $output/01_fastp/log/Phe-{1}-Rep1_{2}_fastp.log "2>&1" ::: 2B 2M 4B 4M ::: L1 L2 L3
        
        echo "[`date`] fastq Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/02_hisat2/" ]; then
        echo
        echo
        echo "[`date`] hisat2"
        echo '-----------------------------------------------'
        mkdir $output/02_hisat2/
        mkdir $output/02_hisat2/log/
	
	source /home/lfgu/bigdata/hitseq/phe_shoot_nanopore/tool/miniconda3/bin/activate hisat2

#	Library Type        Examples                    Tag
#	fr-unstranded       Standard Illumina           FR
#	fr-firststrand      dUTP, NSR, NNSR             RF
#	fr-secondstrand     Ligation, Standard SOLiD   
	
	$base_bin/parallel --jobs 1 hisat2 --max-intronlen 30000 --rna-strandness RF -p 40 -x $hisat2_index -1 $output/01_fastp/Phe-{1}-Rep1_{2}.R1.fastp.fastq -2 $output/01_fastp/Phe-{1}-Rep1_{2}.R2.fastp.fastq -S $output/02_hisat2/Phe-{1}-Rep1_{2}.sam ">" $output/02_hisat2/log/Phe-{1}-Rep1_{2}.hisat2.log "2>&1" ::: 2B 2M 4B 4M ::: L1 L2 L3
	
	conda deactivate	

        echo "[`date`] sam2bam Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/03_sam2bam/" ]; then
        echo
        echo
        echo "[`date`] sam2bam"
        echo '-----------------------------------------------'
        mkdir $output/03_sam2bam/
        mkdir $output/03_sam2bam/log/

	$base_bin/parallel $samtools/samtools view -Sb -o $output/03_sam2bam/Phe-{1}-Rep1_{2}.bam $output/02_hisat2/Phe-{1}-Rep1_{2}.sam ">" $output/03_sam2bam/log/Phe-{1}-Rep1_{2}_sam2bam.log "2>&1" ::: 2B 2M 4B 4M ::: L1 L2 L3
	$base_bin/parallel $samtools/samtools sort $output/03_sam2bam/Phe-{1}-Rep1_{2}.bam -o $output/03_sam2bam/Phe-{1}-Rep1_{2}.sorted.bam ">" $output/03_sam2bam/log/Phe-{1}-Rep1_{2}_bam_sort.log "2>&1" ::: 2B 2M 4B 4M ::: L1 L2 L3
	$base_bin/parallel $samtools/samtools index $output/03_sam2bam/Phe-{1}-Rep1_{2}.sorted.bam ">" $output/03_sam2bam/log/Phe-{1}-Rep1_{2}_bam_index.log "2>&1" ::: 2B 2M 4B 4M ::: L1 L2 L3

	echo "[`date`] sam2bam Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/04_cleanbam/" ]; then
        echo
        echo
        echo "[`date`] cleanbam"
        echo '-----------------------------------------------'
        mkdir $output/04_cleanbam/
        mkdir $output/04_cleanbam/log/
	
	#we can filter the BAM to keep only uniquely mapping reads (aligned concordantly exactly 1 time paired reads).
	$base_bin/parallel $samtools/sambamba view -h -t 3 -f bam -F '"[NH] == 1 and proper_pair"' -o $output/04_cleanbam/Phe-{1}-Rep1_{2}.clean.sorted.bam $output/03_sam2bam/Phe-{1}-Rep1_{2}.sorted.bam ">" $output/04_cleanbam/log/Phe-{1}-Rep1_{2}.cleanbam.log "2>&1" ::: 2B 2M 4B 4M ::: L1 L2 L3
	$base_bin/parallel $samtools/samtools index $output/04_cleanbam/Phe-{1}-Rep1_{2}.clean.sorted.bam ">" $output/04_cleanbam/log/Phe-{1}-Rep1_{2}.index.log "2>&1" ::: 2B 2M 4B 4M ::: L1 L2 L3
	
	echo "[`date`] cleanbam Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/05_TPMCalculator/" ]; then
        echo
        echo
        echo "[`date`] TPMCalculator"
        echo '-----------------------------------------------'
        mkdir $output/05_TPMCalculator/
        mkdir $output/05_TPMCalculator/log/
	
	source /home/lfgu/bigdata/hitseq/phe_shoot_nanopore/tool/miniconda3/bin/activate tpmcalculator

	cd $output/05_TPMCalculator/

	$base_bin/parallel TPMCalculator -g $bamboo_gtf -c 150 -p -b $output/04_cleanbam/Phe-{1}-Rep1_{2}.clean.sorted.bam ">" $output/05_TPMCalculator/log/Phe-{1}-Rep1_{2}.TPMCalculator.log "2>&1" ::: 2B 2M 4B 4M ::: L1 L2 L3
	conda deactivate

	echo "[`date`] TPMCalculator Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/06_featurecounts/" ]; then
        echo
        echo
        echo "[`date`] featurecounts"
        echo '-----------------------------------------------'
	mkdir $output/06_featurecounts/
        mkdir $output/06_featurecounts/log/
	
	$base_bin/parallel $featurecounts/featureCounts -a $bamboo_gtf -p -t exon -g gene_id -o $output/06_featurecounts/Phe-{1}-Rep1_{2}.clean $output/04_cleanbam/Phe-{1}-Rep1_{2}.clean.sorted.bam ">" $output/06_featurecounts/log/Phe-{1}-Rep1_{2}.clean.featurecounts.log "2>&1" ::: 2B 2M 4B 4M ::: L1 L2 L3
	
	echo "[`date`] featureCounts Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/07_deseq2/" ]; then
        echo
        echo
        echo "[`date`] deseq2"
        echo '-----------------------------------------------'
	mkdir $output/07_deseq2/
	mkdir $output/07_deseq2/count/
	mkdir $output/07_deseq2/log/
	
	source /home/lfgu/bigdata/hitseq/phe_shoot_nanopore/tool/miniconda3/bin/activate python3.7
	
	sample1s=(2B 2B 2M 4B)
        sample2s=(2M 4B 4M 4M)
        for ((i=0; i<4; i++));
        do
		python $susupplementary/summarize_count.py --extract_file $susupplementary/summ_Phe_${sample1s[$i]}_vs_${sample2s[$i]}.txt --output $output/07_deseq2/count/count_Phe_${sample1s[$i]}_vs_${sample2s[$i]}.csv > $output/07_deseq2/log/Phe_${sample1s[$i]}_vs_${sample2s[$i]}_count.log 2>&1
	done
	
	conda deactivate
	
	source /home/lfgu/bigdata/hitseq/phe_shoot_nanopore/tool/miniconda3/bin/activate deseq2

	for ((i=0; i<4; i++));
        do
	Rscript $susupplementary/DEseq2.R $output/07_deseq2/count/count_Phe_${sample1s[$i]}_vs_${sample2s[$i]}.csv $output/07_deseq2/deseq2_${sample1s[$i]}_vs_${sample2s[$i]}.csv > $output/07_deseq2/log/Phe_${sample1s[$i]}_vs_${sample2s[$i]}_deseq2.log 2>&1 &
	wait
	done
	
	conda deactivate
	echo "[`date`] DEseq2 Run complete"
        echo '-----------------------------------------------'
fi

#if [ ! -d "$output/08_stringtie/" ]; then
#        echo
#        echo
#        echo "[`date`] stringtie"
#        echo '-----------------------------------------------'
#	mkdir $output/08_stringtie/
#	mkdir $output/08_stringtie/log/
#	
#	$base_bin/parallel $tools/stringtie-2.1.4.Linux_x86_64/stringtie $output/04_cleanbam/Phe-{1}-Rep1_{2}.clean.sorted.bam -G $bamboo_gff --rf -p 10 -o $output/08_stringtie/Phe-{1}-Rep1_{2}.stringtie.gtf ">" $output/08_stringtie/log/Phe-{1}-Rep1_{2}.stringtie.log "2>&1" ::: 2B 2M 4B 4M ::: L1 L2 L3
#	
#	ls $output/08_stringtie/*.gtf > $output/08_stringtie/all_gtf_list.txt
#	$tools/stringtie-2.1.4.Linux_x86_64/stringtie --merge -o $output/08_stringtie/phe_illumina_merge_all.gtf -G $bamboo_gff -p 40 $output/08_stringtie/all_gtf_list.txt > $output/08_stringtie/log/merge.log 2>&1
#
#	echo "[`date`] stringtie Run complete"
#        echo '-----------------------------------------------'
#fi

if [ ! -d "$output/09_rmats/" ]; then
        echo
        echo
        echo "[`date`] rmats"
        echo '-----------------------------------------------'
        mkdir $output/09_rmats/
        mkdir $output/09_rmats/log/
	
	source /home/lfgu/bigdata/hitseq/phe_shoot_nanopore/tool/miniconda3/bin/activate rmats
	
	sample1s=(2B 2B 2M 4B)
        sample2s=(2M 4B 4M 4M)
        for ((i=0; i<4; i++));
        do
		rmats.py --b1 $susupplementary/rmats_Phe_${sample1s[$i]}_rep1.txt --b2 $susupplementary/rmats_Phe_${sample2s[$i]}_rep1.txt --gtf $susupplementary/phe_nanopore_merge.gtf -t paired --readLength 150 --nthread 10 --libType fr-firststrand --novelSS --od $output/09_rmats/Phe_${sample1s[$i]}_vs_${sample2s[$i]} --tmp $output/09_rmats/Phe_${sample1s[$i]}_vs_${sample2s[$i]} > $output/09_rmats/log/Phe_${sample1s[$i]}_vs_${sample2s[$i]}_rmats.log 2>&1 &
	done
	wait
	echo "[`date`] rmats Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/11_bam2bigwig/" ]; then
        echo
        echo
        echo "[`date`] bam2bigwig"
        echo '-----------------------------------------------'
        mkdir $output/11_bam2bigwig/
        mkdir $output/11_bam2bigwig/log/

	source /home/lfgu/bigdata/hitseq/phe_shoot_nanopore/tool/miniconda3/bin/activate deeptools
	
	$base_bin/parallel $samtools/samtools merge --threads 4 $output/11_bam2bigwig/Phe-{}-Rep1.merge.clean.sorted.bam $output/04_cleanbam/Phe-{}-Rep1_L1.clean.sorted.bam $output/04_cleanbam/Phe-{}-Rep1_L2.clean.sorted.bam $output/04_cleanbam/Phe-{}-Rep1_L3.clean.sorted.bam ">" $output/11_bam2bigwig/log/Phe-{}.merge.log "2>&1" ::: 2B 2M 4B 4M
	$base_bin/parallel $samtools/samtools index $output/11_bam2bigwig/Phe-{}-Rep1.merge.clean.sorted.bam ::: 2B 2M 4B 4M

	$base_bin/parallel bamCoverage -b $output/11_bam2bigwig/Phe-{}-Rep1.merge.clean.sorted.bam --normalizeUsing BPM --binSize 5 -o $output/11_bam2bigwig/Phe-{}-Rep1.merge.clean.bw ">" $output/11_bam2bigwig/log/Phe-{}.bamCoverage.log "2>&1" ::: 2B 2M 4B 4M

	echo "[`date`] bam2bigwig Run complete"
        echo '-----------------------------------------------'
fi

