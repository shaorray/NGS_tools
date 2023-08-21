#!/bin/bash

# MIT License

# Copyright (c) [2023] [Rui Shao]

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

echo ""
date
echo "-------------------------------------------------"

# step 0: arguments -----------------------------------

while getopts g:f:n:1:2:b:c:p:m:o:h opts
do case "$opts" in
	g) genome="$OPTARG";;
	f) FOLDER="$OPTARG";;
	n) NAME="$OPTARG";; 
	1) read_1="$OPTARG";;
	2) read_2="$OPTARG";;
	b) barcode1="$OPTARG";;
	c) barcode2="$OPTARG";;
	p) n_core="$OPTARG";;
	m) memory="$OPTARG";;
	o) out_path="$OPTARG";;
	h) echo "
A pipeline for ATAC-seq data processing and quality-check.

Usage ATAC_wrapper.sh -g <hg19/hg38/mm9/mm10> -1 R1.fq.gz -2 R2.fq.gz -o <output path>

Options:	-g	reference genome.
		-f	full path sample folder.
		-n	sample name.
		-1	fastq read 1.
		-2	fastq read 2.
		-b	barcode table 1, in spatial order.
		-c	barcode table 2, in spatial order.
		-p	number of cpu threads.
		-m	allocate memory gb.
		-o	output folder of reports.
"
exit;;
esac
done

NAME=$( echo $NAME | sed 's/\W/_/g' )

echo "genome: "$genome
echo "folder: "$FOLDER
echo "sample: "$NAME
echo "read_1: "$read_1
echo "read_2: "$read_2
echo "barcod1: "$barcode1
echo "barcod2: "$barcode2
echo "n_core: "$n_core
echo "memory: "$memory
echo "report: "$out_path 

cd $FOLDER

mkdir $FOLDER/bam
mkdir $FOLDER/fastq_val
mkdir $FOLDER/log
mkdir $FOLDER/log/fastqc
mkdir $FOLDER/bw
mkdir $FOLDER/peak
mkdir $FOLDER/tmp_data
mkdir $FOLDER/tmp_data/qc_raw_data
mkdir $FOLDER/picard

bowtie2_index="/p300s/liujiang_group/shaor/genomeDir/bowtie_index/mm10"

# step 0: trim barcode --------------------------------
echo ""
echo "-------------------------------------------------"
echo ">>> Trimming barcodes"


if [ ! -f tmp_data/bc_index* ]
then
	# read2
	bbduk.sh \
	in1=$read_1 \
	in2=$read_2 \
	outm1=tmp_data/qc_raw_data/${NAME}_1.fq.gz \
	outm2=tmp_data/qc_raw_data/${NAME}_2.fq.gz \
	k=15 mm=f rcomp=f restrictleft=30 skipr1=t \
	hdist=3 \
	stats=tmp_data/qc_raw_data/${NAME}_stats.read2.txt \
	threads=$n_core \
	literal=CATCGGCGTACGACT

	# linker2
	bbduk.sh \
	in1=tmp_data/qc_raw_data/${NAME}_1.fq.gz \
	in2=tmp_data/qc_raw_data/${NAME}_2.fq.gz \
	outm1=tmp_data/qc_raw_data/${NAME}_linker2_1.fq.gz \
	outm2=tmp_data/qc_raw_data/${NAME}_linker2_2.fq.gz \
	k=30 mm=f rcomp=f restrictleft=60 skipr1=t \
	hdist=3 \
	stats=tmp_data/qc_raw_data/${NAME}_stats.linker2.txt \
	threads=$n_core \
	literal=ATCCACGTGCTTGAGCGCGCTGCATACTTG

	# linker1
	bbduk.sh \
	in1=tmp_data/qc_raw_data/${NAME}_linker2_1.fq.gz \
	in2=tmp_data/qc_raw_data/${NAME}_linker2_2.fq.gz \
	outm1=tmp_data/qc_raw_data/${NAME}_qc_1.fq.gz \
	outm2=tmp_data/qc_raw_data/${NAME}_qc_2.fq.gz \
	k=30 mm=f rcomp=f restrictleft=100 skipr1=t \
	hdist=3 \
	stats=tmp_data/qc_raw_data/${NAME}_stats.linker1.txt \
	threads=$n_core \
	literal=CCCATGATCGTCCGATGCAGTCGTGCCATG

	wait

	/home/shaor/.conda/envs/run/bin/python /p300s/liujiang_group/shaor/scripts/BC_fastq_rename.py \
		-i1 tmp_data/qc_raw_data/${NAME}_qc_1.fq.gz -i2 tmp_data/qc_raw_data/${NAME}_qc_2.fq.gz \
		-b1 $barcode1 -b2 $barcode2 -o1 tmp_data/${NAME}_S1_L001_R1_001.fastq -o2 tmp_data/${NAME}_S1_L001_R2_001.fastq -o3 tmp_data/${NAME}_S1_L001_R3_001.fastq

	wait
	gzip tmp_data/*fastq
fi

wait

# step 1: fastqc --------------------------------------
echo ""
echo "-------------------------------------------------"
echo ">>> fastqc"

if [ ! -f log/fastqc/*html ]
then
	fastqc tmp_data/${NAME}_S1_L001_R1_001.fastq.gz -o log/fastqc
	fastqc tmp_data/${NAME}_S1_L001_R3_001.fastq.gz -o log/fastqc
fi

wait

# step 2a: CellRanger ---------------------------------

echo ""
echo "-------------------------------------------------"
#echo ">>> CellRanger"

if [ ! -d $NAME/outs ]
then
	/p300s/liujiang_group/shaor/opt/cellranger-atac-2.1.0/cellranger-atac count \
		--id=$NAME \
		--reference="/p300s/liujiang_group/shaor/opt/refdata-cellranger-arc-mm10-2020-A-2.0.0" \
	    	--fastqs=tmp_data \
	    	--sample=$NAME \
	    	--localcores=$n_core \
	    	--localmem=$memory

fi

wait

# step 2b: alignment -----------------------------------

echo ""
echo "-------------------------------------------------"
echo ">>> Trimming adapter"

if [ ! -f fastq_val/*gz ]
then	# remove 3' adapters
	trim_galore -j $n_core --paired --quality 20 --length 20 --gzip -o fastq_val tmp_data/${NAME}_S1_L001_R1_001.fastq.gz tmp_data/${NAME}_S1_L001_R3_001.fastq.gz
fi	

step_2_align () {
echo ">>> bowtie2 aligning"

	bowtie2 -x $bowtie2_index -p $n_core --local \
		-1 fastq_val/${NAME}_S1_L001_R1_001_val_1.fq.gz \
		-2 fastq_val/${NAME}_S1_L001_R3_001_val_2.fq.gz | \
		samtools view -bS - | \
		samtools sort - -O "bam" -o bam/${NAME}.bam
	
	samtools index bam/${NAME}.bam
	
	samtools flagstat bam/${NAME}.bam > log/${NAME}.flagstat.log
	samtools idxstats bam/${NAME}.bam > log/${NAME}.idxstats.log
}

if [ ! -f bam/*bai ]
then
	step_2_align
fi

wait

# step 3: peak calling --------------------------------
echo ""
echo "-------------------------------------------------"
echo ">>> peak calling"

if [ ! -f peak/*narrowPeak ]
then
	macs2 callpeak -t bam/${NAME}.bam -f BAMPE -g mm -n ${NAME}"_ATAC" -B -q 0.05 \
 --nomodel --keep-dup all --shift -100 --extsize 200 --outdir peak \
 2> log/${NAME}.ATAC_nm_MACS2.log
	
	wait
	
	samtools sort -n bam/${NAME}.bam > bam/${NAME}_sorted.bam
	wait
	
	Genrich -t bam/${NAME}_sorted.bam -o peak/${NAME}.genrich.narrowPeak -p 0.03 -j -y -r -e chrM -v
	
	wait

	awk '{print $1"\tGenrich\tpeak\t"($2+1)"\t"$3"\t"$5"\t"$6"\t.\tname=" $4";signalValue=" $7";pValue=" $8";qValue=-1;peak=" $10}' peak/${NAME}.genrich.narrowPeak > peak/${NAME}.genrich.narrowPeak.gtf

fi

# step 4: reads coverage ------------------------------
echo ""
echo "-------------------------------------------------"
echo ">>> bigwig conversion"

if [ ! -f bw/*bw ]
then
	bamCoverage -b bam/${NAME}.bam -o bw/${NAME}.bw --binSize 10 -p $n_core
fi

wait 

# step 5: insertion size distrubution -----------------
echo ""
echo "-------------------------------------------------"
echo ">>> insertion size"


if [ ! -f picard/* ]
then
	picard CollectInsertSizeMetrics \
		I=bam/${NAME}.bam O=picard/$i.insert_size_metrics.txt \
		H=picard/${NAME}.insert_size_histogram._raw.pdf M=0.5
fi

wait

# step 6: reads quality and distribution --------------
echo ""
echo "-------------------------------------------------"
echo ">>> R quality control"

#/home/shaor/.conda/envs/run/bin/Rscript /p300s/liujiang_group/shaor/scripts/ATAC_QC.R -b bam/${NAME}.bam -p peak/${NAME}".genrich.narrowPeak" -g $genome -d $out_path

# step 7: get pixel read counts ------------------------
echo ""
echo "-------------------------------------------------"
echo ">>> run HTSeq count"

if [ ! -f bam/${NAME}_tagged.bam ]
then
	/home/shaor/.conda/envs/run/bin/python /p300s/liujiang_group/shaor/scripts/BC_add_bam_tag.py -i bam/${NAME}.bam  -o bam/${NAME}_tagged.bam
fi

wait

# add barcodes to BAM header as RG tag (readGroup)
if [ ! -f feature* ]
then
	samtools view -H bam/${NAME}_tagged.bam > bam/header.sam
	bc_n=$(wc -l < $barcode1)
	/home/shaor/.conda/envs/run/bin/python /p300s/liujiang_group/shaor/scripts/BC_sam_RG_header.py -i $bc_n -o bam/bc_RG
	cat bam/bc_RG bam/header.sam > bam/header_bc.sam
	samtools reheader bam/header_bc.sam bam/${NAME}_tagged.bam > bam/${NAME}_tagged_bc.bam
	samtools index bam/${NAME}_tagged_bc.bam
	
	samtools sort bam/${NAME}_tagged.bam > bam/${NAME}_sort_tagged.bam
	samtools index bam/${NAME}_sort_tagged.bam
	
	wait
	
	featureCounts -p --countReadPairs --primary -Q 5 -O -T $n_core --byReadGroup \
	 -a peak/${NAME}.genrich.narrowPeak.gtf -F GFF -t peak -g name \
	 -o featureCounts_output.txt bam/${NAME}_tagged_bc.bam
fi
