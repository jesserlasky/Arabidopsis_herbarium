## For herbaria (ancient) samples

## LeeHom script _ LeeConcat.sh

#!/bin/sh

#  LeeConcat.sh
#
#
# setting variables
SAMPLES=$(cat SamplesNameFile.txt)

RAWREADS=$(echo "./PATH_TO_SAMPLES")
LEEHDIR=$(echo "./LeeHom")
READSDIR=$(echo "./PATH_TO_READS")
#
#
# adapter trimming and merge reads
for i in $SAMPLES
do
    	echo "LeeHom sample $i starting"
    	/usr/local/bin/leeHom/src/leeHom --ancientdna -fq1 ${RAWREADS}/${i}_R1_001.fastq.gz -fq2 ${RAWREADS}/${i}_R2_001.fastq.gz -fqo ${LEEHDIR}/${i}
    	echo "LeeHom sample $i done"
done &> >(tee ${READSDIR}/file.out)
#
#
##### concatenate files were reads are unmerged but properly paired with this file R1R2+R1+R2
for i in $SAMPLES
do
    	echo "cat sample $i"    
    	cat ${LEEHDIR}/${i}.fq.gz ${LEEHDIR}/${i}_r1.fq.gz ${LEEHDIR}/${i}_r2.fq.gz> ${READSDIR}/${i}_all.fq.gz
    	echo "cat done for sample $i"
done
#

## Mapping to the genome_mapping.sh
#!/bin/sh

#  mapping.sh
#  
#
# setting variables
ALL_SAMPLES=$(cat SamplesNameFile.txt)
REFDIR=$(echo "./PATH_TO_GENOME/ ")
READSDIR=$(echo "./PATH_TO_READS")
ALNDIR=$(echo "./PATH_TO_ALN_OUTPUT")
MAINDIR=$(echo "./PATH_TO_")
FINALDIR=$(echo "./FilterMapALL")
#
#
# indexing, dictionary
bwa index ${REFDIR}/Athaliana_TAIR10.fa
java -jar /usr/bin/picard-2.8.2/picard-2.8.2.jar CreateSequenceDictionary REFERENCE=${REFDIR}/Athaliana_TAIR10.fa OUTPUT=${REFDIR}/At
haliana_TAIR10.dict
samtools faidx Athaliana/Araport11/assembly/Athaliana_TAIR10.fa
#
# mapping aln (-n, -o, -l values are specific for aDNA) -t is number of cores
for i in $ALL_SAMPLES
do
	bwa aln -n 0.01 -o 2 -l 16500 -t 7 ${REFDIR}/Athaliana_nc_TAIR10.fa ${READSDIR}/${i}_all.fq.gz > ${ALNDIR}/${i}.sai
done &> >(tee ${MAINDIR}/alnnor.out)
#
## mapping samse
for i in $ALL_SAMPLES
do
	echo "Starting bwa mem for sample $i ..."
	RG=$(echo "@RG\tID:${i}\tSM:${i}\tPL:Illumina\tLB:${i}")
	bwa samse -r ${RG} ${REFDIR}/Athaliana_nc_TAIR10.fa ${ALNDIR}/${i}.sai ${READSDIR}/${i}_all.fq.gz > ${MAINDIR}/${i}.sam
	echo "samse $i done"
done &> >(tee ${MAINDIR}/samsenor.out)
#
#sam to bam
for i in $ALL_SAMPLES
do
	samtools view -bhS ${MAINDIR}/${i}.sam > ${MAINDIR}/${i}.bam

	rm ${MAINDIR}/${i}.sam

	echo "'sam to bam' and 'remove sam' for sample $i done!"

done &> >(tee ${MAINDIR}/samtobamnor.out)
#
# sort
for i in $ALL_SAMPLES
do
	echo "$i sorting"

	samtools sort -T ${MAINDIR}/${i}.sorted -o ${MAINDIR}/${i}.sorted.bam ${MAINDIR}/${i}.bam

	echo "$i sorted"
done &> >(tee ${MAINDIR}/sortnor.out)
#
# checking mapping stats
for i in $ALL_SAMPLES
do
	samtools flagstat ${MAINDIR}/${i}.sorted.bam
done &> >(tee ${MAINDIR}/mappingstatsnor.out)
#
# remove duplicates
for i in $ALL_SAMPLES
do
	echo "$i removing duplicates"
	samtools rmdup -S ${MAINDIR}/${i}.sorted.bam ${MAINDIR}/${i}_dedup.bam
	echo "$i duplicates removed"
done &> >(tee ${MAINDIR}/dedupnor.out)
#
# check mapping stats
for i in $ALL_SAMPLES
do
	samtools flagstat ${MAINDIR}/${i}_dedup.bam
done &> >(tee ${MAINDIR}/statsmapdedupnor.out)
#  filter quality-length (min size 30bp and quality 30) and index bam files (IMPORTANT that the bam files are indexed or GATK won't work)
for i in $ALL_SAMPLES
do
	samtools view -b -m 30 -q 30 ${MAINDIR}/${i}_dedup.bam -o ${FINALDIR}/${i}.bam
	samtools index ${FINALDIR}/${i}.bam
	echo "$i FQ filtered and indexed"
done &> >(tee ${FINALDIR}/HQnor.out)

