#!/bin/bash
set -e

#REQUIRED TOOLS: BWA, SAMTOOLS, PICARD, GATK, BCFTOOLS, vcf2consensus.pl, Muscle3.8.
#V1 Basic BWA alignment to given reference and consensus generation (edited from Richard Ellis IterMap.sh)
#V2 Added indel realignment using GATK and read group renaming with Picard
#V3 Removed need to use text entered to name samples, changed so sample name taken from read and ref files
#V4 Added counter to automatically run iterations, added iteration requirement to command
#V5 Added LogRun command and use of output directory
#V6 Added seperation of consensus calling to remove duplicates from final consensus calling bam
#V7 Removed A switch from bcftools call command as was forcing all SNPs into consensus no matter how prevalent
#B-version Removed sym-link so need to be in same directory as ref and reads. Removed output directory options
#VB8 added multi-threading option
#VB9 Re-added A switch to bcftools during cycling to prevent N incorperation
#VB10 Added default and command options for threads and iterations
#VB11 COMMENTED OUT SECTIONS DELETED, full version archived and saved, version A also archived. Also GATK and picard stuuf removed (file names altered accoridngly.


# set defaults for the options

iter=1
threads=$(grep -c ^processor /proc/cpuinfo)
count=1
Start=$(date +%s)

# parse the options
while getopts 'i:t:' opt ; do
  case $opt in
    i) iter=$OPTARG ;;
    t) threads=$OPTARG ;;
  esac

done

# skip over the processed options
shift $((OPTIND-1))

# check for mandatory positional parameters and define what position is what
if [ $# -lt 3 ]; then
	echo "Usage: $0 [-i desired iterations] [-t threads] {path to reference.fa} {path to read1} {path to read2}"

	exit 1
fi

REF="$1"
READ1="$2"
READ2="$3"

# define names of sample and ref removing the filoe extensions from their names
sfile1=$(basename "$READ1")
sfile2=$(basename "$READ2")
samplename=${sfile1%%_*}

ref=$(basename "$REF")
refname=${ref%%.*}

#Start counting and iteration cycling
while (($count <= $iter))
do

#Create logrun file and create the LogRun command
echo -e "$now\n\tbulkRGA_$ref.\n the following commands were run:\n" > refguidealign_"$refname".log
LogRun(){
echo -e "\n$@" >> refguidealign_"$refname".log
eval "$@"
}

#index the reference sequence
LogRun bwa index "$REF"
LogRun samtools faidx "$REF"

#create reference dictionary
#if [ -f "${REF%%.*}.dict" ]
#then echo "Using existing reference dictionary"
#else echo "Creating reference dictionary"
#	LogRun java -jar /home/p991369/bin/picard.jar CreateSequenceDictionary R="$REF" O=${REF%.*}.dict
#fi

#align to reference and convert sam to bam, then sort
LogRun bwa mem -t "$threads" -T 10 -k 18 -B 4 -O 6 "$REF" "$READ1" "$READ2" | samtools sort -@ "$threads" -O bam > "$samplename"-"$refname"-ITER"$count"sorted.bam

#assign read groups otherwise GATK gets arsy
#LogRun java -jar /home/p991369/bin/picard.jar AddOrReplaceReadGroups I="$samplename"-"$refname"-ITER"$count"sorted.bam O="$samplename"-"$refname"-ITER"$count"assign.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

#index the bam
LogRun samtools index "$samplename"-"$refname"-ITER"$count"sorted.bam

#realign consensus to sort out indels
#LogRun java -jar /opt/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt "$threads" -R "$REF" -I "$samplename"-"$refname"-ITER"$count"assign.bam -o indel"$samplename"-"$refname"-ITER"$count".list

#LogRun java -jar /opt/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R "$REF" -I "$samplename"-"$refname"-ITER"$count"assign.bam -targetIntervals indel"$samplename"-"$refname"-ITER"$count".list -maxReads 50000 -o "$samplename"-"$refname"-ITER"$count"realigned.bam

#if on the last iteration duplicates are marked and removed to allow the creation of a 'clean' bam from which to create a consensus

if [ $count == $iter ]; then

#LogRun java -jar /home/p991369/bin/picard.jar MarkDuplicates INPUT="$samplename"-"$refname"-ITER"$count"realigned.bam OUTPUT="$samplename"-"$refname"-ITER"$count"cleaned.bam REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics"$count".txt

#create pileup file then vcf used to make consensus
LogRun samtools mpileup -L 10000 -Q 0 -AEupf "$REF" "$samplename"-"$refname"-ITER"$count"sorted.bam | bcftools call -c  > "$samplename"-"$refname"-ITER"$count".vcf

#make consensus and change header of the fasta file
LogRun vcf2consensus.pl consensus -d 2 -Q 0 -f "$REF" "$samplename"-"$refname"-ITER"$count".vcf | sed '/^>/ s/-ITER[0-9]//;/^>/ s/$/'-ITER"$count"'/' > "$samplename"-"$refname"-ITER"$count"-consensus.fa

#create mapping stats
LogRun samtools flagstat "$samplename"-"$refname"-ITER"$count"sorted.bam > "$samplename"-"$refname"-ITER"$count"stats.txt

LogRun samtools index "$samplename"-"$refname"-ITER"$count"sorted.bam

else

#create pileup file then vcf used to make consensus
LogRun samtools mpileup -L 10000 -Q 0 -AEupf "$REF" "$samplename"-"$refname"-ITER"$count"sorted.bam | bcftools call -Ac > "$samplename"-"$refname"-ITER"$count".vcf

#make consensus and change header of the fasta file
LogRun vcf2consensus.pl consensus -d 2 -Q 0 -f "$REF" "$samplename"-"$refname"-ITER"$count".vcf | sed '/^>/ s/-ITER[0-9]//;/^>/ s/$/'-ITER"$count"'/' > "$samplename"-"$refname"-ITER"$count"-consensus.fa

#create mapping stats
LogRun samtools flagstat "$samplename"-"$refname"-ITER"$count"sorted.bam > "$samplename"-"$refname"-ITER"$count"stats.txt

fi

#Assign new reference sequence
REF="$samplename"-"$refname"-ITER"$count"-consensus.fa

#Add to the iteration count/ end cycling
((count=count+1))

done

#remove un-needed files
#rm *assign.bam
#rm *sorted.bam
#rm *assign.bam.bai
#rm *.amb
#rm *.ann
#rm *.bwt
#rm *.fai
#rm *.pac
#rm *.sa
#rm unaligned.fa
#rm *.list
#rm *.dict

End=$(date +%s)
TimeTaken=$((End-Start))
echo  | awk -v D=$TimeTaken '{printf "Performed '$iter' mapping iterations in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'





