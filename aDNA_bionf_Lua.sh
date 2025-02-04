## For herbaria (ancient) samples

## LeeHom script name  LeeConcat.sh

#!/bin/sh

#  LeeConcat.sh
# 
# Created by Lua Lopez 
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
    	echo "cat done for sample $i”
done
#

## Mapping to the genome script name mapping.sh
#!/bin/sh

#  mapping.sh
#  
# Created by Lua Lopez 
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


## For non-herbaria (modern) samples

### 1001G and Durvasula et al. 2017 samples were downloaded from NCBI and the European Nucleotide Archive/Sequence Read Archive

### Trimming script name trimmomatic.sh

#### Trimmomatic
#!/bin/sh
#  trimmomatic.sh
#  
#
#  Created by Lua Lopez  
#
# setting variables
ALL_SAMPLES=$(cat SamplesNameFile.txt)
RAWREADS=$(echo "./PATH_TO_RAW_READS/")
READSDIR=$(echo "./PATH_TO_TRIMMED_READS/")
#
for i in $ALL_SAMPLES
do
echo “Trimmomatic for sample $i starting”
java -jar /home/lxl326/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 ${RAWREADS}/${i}_1.fastq.gz ${RAWREADS}/${i}_2.fastq.gz ${READSDIR}/${i}_1_paired.fq.gz ${READSDIR}/${i}_1_unpaired.fq.gz ${READSDIR}/${i}_2_paired.fq.gz ${READSDIR}/${i}_2_unpaired.fq.gz ILLUMINACLIP:/home/lxl326/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:3:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:100
echo "Trimmomatic done for sample $i"
done &> >(tee ./Trimmomatic.out)


## mapping script for non-herbaria samples mapping.sh
#!/bin/sh
#  mapping.sh
#  
#
#  Created by Lua Lopez  
#
# setting variables
ALL_SAMPLES=$(cat SamplesNameFile.txt)
REFDIR=$(echo "./PATH_TO_GENOME/")
READSDIR=$(echo ""./PATH_TO_TRIMMED_READS/")
BWADIR=$(echo "./PATH_TO_BWA_OUTPUT/")
FINALDIR=$(echo "./PATH_TO_FILTERED_READS")
#
#
# indexing, dictionary
bwa index ${REFDIR}/Athaliana_TAIR10.fa
java -jar /usr/bin/picard-2.8.2/picard-2.8.2.jar CreateSequenceDictionary REFERENCE=${REFDIR}/Athaliana_TAIR10.fa OUTPUT=${REFDIR}/Athaliana_TAIR10.dict
samtools faidx /media/elgon/Lua/aDNA_2020/Athaliana/Araport11/Athaliana_TAIR10.fa
#
## mapping mem
for i in $ALL_SAMPLES
do
	echo "Starting bwa mem for sample $i ..."
	RG=$(echo "@RG\tID:${i}\tSM:${i}\tPL:Illumina\tLB:${i}")
	bwa mem -t 17 -R ${RG} ${REFDIR}/Athaliana_nc_TAIR10.fa ${READSDIR}/${i}_1_paired.fq.gz ${READSDIR}/${i}_2_paired.fq.gz > ${BWADIR}/${i}.sam
done &> >(tee ${BWADIR}/mem.out)
#
#sam to bam
for i in $ALL_SAMPLES
do
	echo “ starting ‘sam to bam' and 'remove sam' for sample $i”
samtools view -bhS ${BWADIR}/${i}.sam > ${BWADIR}/${i}.bam
	rm ${BWADIR}/${i}.sam
	echo "'sam to bam' and 'remove sam' for sample $i done!"
done &> >(tee ${SAMDIR}/samtobam.out)
#
# sort
for i in $ALL_SAMPLES
do
	echo "$i sorting"
	samtools sort -T ${BWADIR}/${i}.sorted -o ${BWADIR}/${i}.sorted.bam ${BWADIR}/${i}.bam
	echo "$i sorted"
done &> >(tee ${BWADIR}/sort.out)
#
# checking mapping stats
for i in $ALL_SAMPLES
do
    	echo "$i checking mapping stats"
	samtools flagstat ${BWADIR}/${i}.sorted.bam
	echo "$i done mapping stats"
done &> >(tee ${BWADIR}/mappingstats.out)
#
# remove duplicates
for i in $ALL_SAMPLES
do
	echo "$i removing duplicates"
	samtools rmdup -S ${BWADIR}/${i}.sorted.bam ${BWADIR}/${i}_dedup.bam
	echo "$i duplicates removed"
done &> >(tee ${BWADIR}/dedup.out)
#
# check mapping stats
for i in $ALL_SAMPLES
do
    	echo "$i starting dedup stats"
	samtools flagstat ${BWADIR}/${i}_dedup.bam
	echo "$i done dedup stats"
done &> >(tee ${BWADIR}/statsmapdedup.out)
#
#  filter quality-length (min size 30bp and quality 30) and index bam files (IMPORTANT that the bam files are indexed or GATK won't work)
for i in $ALL_SAMPLES
do
	Echo “starting $i quality filter and index”
samtools view -b -m 30 -q 30 ${BWADIR}/${i}_dedup.bam -o ${FINALDIR}/${i}.bam
	samtools index ${FINALDIR}/${i}.bam
	echo "$i FQ filtered and indexed"
done &> >(tee ${FINALDIR}/HQ.out)


### SNP calling process for all the samples. 

# Varscan individual genotype calling 

## Varscan script name Varscan.sh NOTE samtools and varscan does not allow you to specify multiple cores. To speed it up the process, the file containing the samples names was divided into SEVERAL files. 

#!/bin/sh
#  Varscan.sh
#  
#  Created by Lua Lopez  
#
# setting variables
ALL_SAMPLES=$(cat SamplesNameFile.txt)
REFDIR=$(echo "./PATH_TO_GENOME/")
MAINDIR=$(echo "./PATH_TO_BAM_FILES")
FINALDIR=$(echo "./PATH_TO_VAR_FILES")
#
#
# sort
for i in $ALL_SAMPLES
do
	echo "$i into varscan"

	samtools mpileup -f ${REFDIR}/Athaliana_TAIR10.fa -q 20 ${MAINDIR}/${i}.bam | java -jar /home/lxl326/VarScan.v2.3.9.jar mpileup2snp --min-coverage 2 --min-reads2 2 > ${FINALDIR}/${i}.vars.txt

	echo "$i donE in varscan"
done &> >(tee ${FINALDIR}/varscan.out)
#

# check the coverage in all the samples. Script Qualtile.pl NOTE: to be run in the folder where all the vars.txt files are. Quantile.pl has to be inside of the folder! 

for i in *.vars.txt ; do cat $i | awk 'NR>1' | cut -f5 | cut -f2 -d: | ./quantile.pl > `basename $i .txt `.quantiles.txt ; done
cat *quantiles.txt | awk 'NF==2' > allSamps.varQuantiles.txt

# where the Quantile.pl is 
#!/usr/bin/perl

#pipe in a column of integers, and float args for all desired points. return quantiles for the floats

#Logan Kistler, Smithsonian Institution

@points = @ARGV;
unless (@points) {
    @points = (.1, .25, .5, .75, .9, .95, .99, .999, 1);
}
@points = sort {$a <=> $b} @points;

while ($line = <STDIN>) {
    chomp $line;
    $hold{$line}++;
    $count++;
}

@keys = sort {$a <=> $b} keys %hold;
foreach $i (@points) {
    $n = int ($count*$i);
    push @set, $n;
}
while (@set>0) {
    $target = shift @set;
    until ($run >= $target) {
   	 $pull = shift @keys;
   	 $run += $hold{$pull};
    }
    $t = shift @points;
    print "$t\t$pull\n";
}

#To create the plot and visualize the results In R:
d <- read.table("allSamps.varQuantiles.txt")
library(lattice)
stripplot(d[,2]~as.character(d[,1]),jitter=T)

### once the threshold has been determined (let’s say 0.9 for the example) we create the bed file 
# #Logan Kistler, Smithsonian Institution, 2018

for i in *.vars.txt ; do f=`basename $i .txt `; q=`awk '$1==0.9 {print $2}' $f.quantiles.txt` ; cat $i | perl -e '$maxcov='$q' ; $line = <> ; while ($line = <>) {chomp $line ; $cv = (split /\s+/, $line)[4] ; $cv = (split /\:/, $cv)[1] ; if ($cv <= $maxcov) {@d = split /\s+/, $line ; $chr = $d[0] ; $pos = $d[1] ; $lopos = $pos-1 ; print "$chr\t$lopos\t$pos\n"}}' ; done | perl -e 'while ($line = <>) {chomp $line ; $hold{$line}++ ; if ($hold{$line}==1) {print $line."\n"}}' | bedtools sort -i - > mySNPs.bed

### mpileup bash script 
#!/bin/sh
# mpileup.sh
#
#
# Created by Lua Lopez
#
# setting variables
SAMPLES=$(cat SamplesNameFile.txt)
REFDIR=$(echo "./PATH_TO_GENOME/")
BEDDIR=$(echo "./PATH_TO_FOLDER_BED_FILE/ ")
READSDIR=$(echo "./PATH_TO_FOLDER_BAM_FILES/")
MPILEDIR=$(echo "./PATH_OUTPUT_FOLDER/")
for i in $SAMPLES
do
    	echo "mpileup for sample $i"    
samtools mpileup -l ${BEDDIR}/mySNPs.bed -aa -f ${REFDIR}/Athaliana_TAIR10.fa ${READSDIR}/${i}.bam > ${MPILEDIR}/${i}.mpileup.txt
    	echo "mpileup done for sample $i"
done &> >(tee ${MPILEDIR}/mpileup.out)

# script name pile.grab NOTE: using max coverage of 70, that is when in the graph points (0.9) start dispersing and min coverage 2

#!/bin/sh
# Pile_Grab.sh
#
#
# Created by Lua Lopez
#
# setting variable
SAMPLES=$(cat SamplesNameFile.txt)
FINALDIR=$(echo "./hap_grab_24")
for i in $SAMPLES
do
    	echo "hapgrab for sample $i"
            	cat ${i}.mpileup.txt | perl stripPile_v2.pl | perl hapGrab.pl /dev/stdin ${i} 2 70 > $FINALDIR/${i}.min2x.max70.hap.txt
    	echo "hapgrab done for sample $i"
done &> >(tee Pile_Grab.out)

### Where script stripPile_v2.pl

#!/usr/bin/perl
#Logan Kistler, Smithsonian Institution, 2018

while ($line = <STDIN>) {
	chomp $line;
	@d = split /\s+/, $line;
	$ping = $d[4];
	$ref = $d[2];
	$ping =~ s/[\.\,]/$ref/g;
	$ping =~ s/[^acgtACGT]//g;
	$ping =~ tr/[acgt]/[ACGT]/;
	$len = length $ping;
	$ping =~ s/A//g;
	$a = $len - (length $ping);
	$len = length $ping;
	$ping =~ s/C//g;
	$c = $len - (length $ping);
	$len = length $ping;
	$ping =~ s/G//g;
	$g = $len - (length $ping);
	$t = length $ping;
	$tot = $a + $c + $g + $t;
	pop @d;
	pop @d;
	pop @d;
	pop @d;
	$set = join "\t", @d;
	print "$set\t$tot\t$a\t$c\t$g\t$t\n";
}

=pod
#!/usr/bin/perl
#
#Takes a bed file of SNP sites (chr, beg, fin) and a pileup file run through stripPile.pl, and adds any missing sites from the bed to make a complete pileup
#

($pile, $bed, $out) = @ARGV;

die "\nUsage: perl $0 infile.pile markers.bed outfile\n\n" unless (@ARGV==3);

open IN, $pile;
while ($line = <IN>) {
	chomp $line;
	($chr, $pos, $n, $a, $c, $g, $t) = split /\s+/, $line;
	$ping = "$chr.$pos";
	$ct = $a + $c + $g + $t;
	$line = join "\t", ($chr, $pos, $ct, $a, $c, $g, $t);
	$hold{$ping} = $line;
}
close IN;

open BED, $bed;
open OUT, ">$out";
while ($line = <BED>) {
	chomp $line;
	($chr, $n, $pos) = split /\s+/, $line;
	$ping = "$chr.$pos";
	$hold{$ping} = "$chr\t$pos\t0\t.\t.\t.\t." unless (exists ($hold{$ping}));
	print OUT $hold{$ping}."\n";
}
## where script hapGrab.pl
#Logan Kistler, Smithsonian Institution, 2018

#!/usr/bin/perl

($in, $id, $cvg, $max) = @ARGV;
$max = 1e12 unless ($max);
open IN, $in;
print "$id\n";
while ($line = <IN>) {
	($chr, $beg, $tot, $a, $c, $g, $t) = split /\s+/, $line;
	@d = split /\s+/, $line;
	if ($a eq "." or $tot > $max) {
		print "0\n";
		next;
	}
	$a = 0 if ($a < $cvg);
	$c = 0 if ($c < $cvg);
	$g = 0 if ($g < $cvg);
	$t = 0 if ($t < $cvg);
	$seq = "";
	$seq .= "A"x$a;
	$seq .= "C"x$c;
	$seq .= "G"x$g;
	$seq .= "T"x$t;
	if ((length $seq) < $cvg) {
		print "0\n";
		next;
	}
	$ping = int(rand(length $seq));
	$ping = substr $seq, $ping, 1;
	print "$ping\n";
}
### creating the pseudohaplotypes using the script newHapToPlink_v2.pl
screen perl newHapToPlink_v2.pl --bed mySNPs.bed --out outfile_stem --biallel *.min2x.max70.hap.txt

## where script newHapToPlink_v2.pl
#Logan Kistler, Smithsonian Institution, 2018
#!/usr/bin/perl
#
#Take a bed file used to make hapGrab.pl files (1), an out stem (2), and a set of hapGrab.pl files (3...)
#make a tped and tfam file.

use Getopt::Long;

GetOptions (
	'bed=s' => \$bed,
	'out=s' => \$out,
	'invar' => \$invar,	#invoke this to filter out invariant sites
	'biallel' => \$biallel,	#enforce max biallelic (does not invoke invar)
	'last=s' => \$last,	#give chr,pos for the last SNP to check. for debugging.
);


#$haps = $ARGV[0];

@files = @ARGV;
#@files = split /\,/, $haps;
#print @files."\n";
#print $hapfiles."\n";
#print $haps."\n";
#print $bed."\n";
#print $invar."\n";
if ($last) {
	($lchr, $lpos) = split /\,/, $last;
}
	

#die "\nUsage: perl $0 snps.bed outstem file1 file2 file3 ...\n\n" unless (@ARGV);

open BED, $bed;
#open OUT, ">$out.sites.tmp";
#while ($line = <IN>) {
#	chomp $line; 
#	@d = split /\s+/, $line;
#	$chr = shift @d; 
#	$null = shift @d;
#	$pos = shift @d;
#	$col2 = "$chr\_$pos";
#	print OUT "$chr $col2 0 $pos\n";
#}
#close OUT;

open TFAM, ">$out.tfam";
open OUT, ">$out.tped";
$set = join " ", @files;
open IN, "paste $set|";
$line = <IN>;
@d = split /\s+/, $line;
foreach $i (@d) {
	print TFAM "$i 1 0 0 0 -9\n";
}
close TFAM;
while ($line = <IN> and $b = <BED>) {
	if ($. % 100000 == 0) {
		print STDERR $.." sites handled\n";
	}
	chomp $b; 
	@d = split /\s+/, $b;
	$chr = shift @d; 
	$null = shift @d;
	$pos = shift @d;
	if ($last) {
		exit if ($lchr eq $chr and $lpos eq $pos);
	}
	$col2 = "$chr\_$pos";
	$newline ="$chr $col2 0 $pos";
	chomp $line;
	@d = split /\s+/, $line;
	if ($invar or $biallel) {
		undef %hold;
	}
	foreach $i (@d) {
		$hold{$i}++;
		$newline .= " $i $i";
#		print OUT " $i $i";
	}
	delete $hold{0};
	if ($invar) {
		next if (keys %hold <= 1);
	}
	if ($biallel) {
		next if (keys %hold > 2);
	}
#	$newline =~ s/^\ //;
	print OUT "$newline\n";	
}
close OUT;

## creating vcf file with no filtering option (those can be added later during each analysis) using plink 

./plink --tfile file_stem --out ../../Final_VCF_24/snp_noFQ --recode vcf


