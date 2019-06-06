while getopts h:o:1:2:e:g:p:m: option
do
case "${option}"
in
h) hictoolPath=${OPTARG};; # The path where are the HiCtool scripts with the final trailing slash.
o) outputPath=${OPTARG};; # The path where to save the output files.
1) fastq1=${OPTARG};; # The fastq file with the first reads of the pairs.
2) fastq2=${OPTARG};; # The fastq file with the second reads of the pairs.
e) restrictionEnzyme=${OPTARG};; # The restriction enzyme or enzymes passed between square brackets (example: [enzyme1,enzyme2]).
g) genomeIndex=${OPTARG};; # The Bowtie2 genome indexes of the reference genome (only filename without extension).
p) threads=${OPTARG};; # The number of parallel threads to use for alignment and pre-processing. The more the fastest the process.
m) max_lines=${OPTARG};; # The maximum number of lines per each temporary fastq file in order to avoid memory errors. Each temporary file is processed by a separate processor if multiple threads are used.
esac
done

echo "Start data preprocessing: $(date)"

checkMakeDirectory(){
if [ ! -e "$1" ]; then
mkdir -p "$1"
fi
}
checkMakeDirectory $outputPath
cd $outputPath

echo -n "Calculating total lines of the fastq files ... "
fastq_lines=`wc -l $fastq1 | awk '{print $1}'`
echo "Done!"

if [ -z $max_lines ]
then
	echo "max_lines not declared."
	python $hictoolPath"HiCtool_pre_truncation.py" -i [$fastq1,$fastq2] -e $restrictionEnzyme -p $threads

	tot_reads_1=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_reads_2=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq2%%.*}_log.txt")

	tot_trunc_reads_1=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_trunc_reads_2=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq2%%.*}_log.txt")

	perc1=$(awk -v n1=$tot_trunc_reads_1 -v ntot1=$tot_reads_1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
	perc2=$(awk -v n1=$tot_trunc_reads_2 -v ntot1=$tot_reads_2 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)

	read_length=$(awk 'FNR == 1 {print $2}' "${fastq1%%.*}_log.txt")

	printf $fastq1"\n"$tot_reads_1" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_1" ("$perc1"%%)\n\n" > pre_truncation_log.txt
	printf $fastq2"\n"$tot_reads_2" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_2" ("$perc2"%%)" >> pre_truncation_log.txt

	rm "${fastq1%%.*}_log.txt"
	rm "${fastq2%%.*}_log.txt"

elif ! [ -z $max_lines ] && [ $max_lines -ge $fastq_lines ]
then
	echo "max_lines not consider because greater that the total lines of the fastq file."
	python $hictoolPath"HiCtool_pre_truncation.py" -i [$fastq1,$fastq2] -e $restrictionEnzyme -p $threads

	tot_reads_1=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_reads_2=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq2%%.*}_log.txt")

	tot_trunc_reads_1=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_trunc_reads_2=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq2%%.*}_log.txt")

	perc1=$(awk -v n1=$tot_trunc_reads_1 -v ntot1=$tot_reads_1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
	perc2=$(awk -v n1=$tot_trunc_reads_2 -v ntot1=$tot_reads_2 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)

	read_length=$(awk 'FNR == 1 {print $2}' "${fastq1%%.*}_log.txt")

	printf $fastq1"\n"$tot_reads_1" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_1" ("$perc1"%%)\n\n" > pre_truncation_log.txt
	printf $fastq2"\n"$tot_reads_2" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_2" ("$perc2"%%)" >> pre_truncation_log.txt

	rm "${fastq1%%.*}_log.txt"
	rm "${fastq2%%.*}_log.txt"

elif ! [ -z $max_lines ] && [ $max_lines -lt $fastq_lines ]
then
	echo -n "Using max_lines to split the fastq files ... "
	if (( $max_lines % 4 )) ; then
		max_lines=`expr $max_lines - $(($max_lines % 4))`
	fi
	# Splitting the first fastq file
	k=$max_lines
	count=1
	while [ $k -lt $fastq_lines ]
	do
		start=`expr $k - $max_lines + 1`
		quit=`expr $k + 1`
		sed -n "$start,"$k"p;"$quit"q" $fastq1 > "${fastq1%%.*}_temp_"$count".fastq"
		count=`expr $count + 1`
		k=`expr $k + $max_lines`
	done
	start=`expr $k - $max_lines + 1`
	sed -n "$start,"$fastq_lines"p" $fastq1 > "${fastq1%%.*}_temp_"$count".fastq"

	# Splitting the second fastq file
	k=$max_lines
	count=1
	while [ $k -lt $fastq_lines ]
	do
		start=`expr $k - $max_lines + 1`
		quit=`expr $k + 1`
		sed -n "$start,"$k"p;"$quit"q" $fastq2 > "${fastq2%%.*}_temp_"$count".fastq"
		count=`expr $count + 1`
		k=`expr $k + $max_lines`
	done
	start=`expr $k - $max_lines + 1`
	sed -n "$start,"$fastq_lines"p" $fastq2 > "${fastq2%%.*}_temp_"$count".fastq"
	echo "Done!"

	# Generate list of temporary files to pass to pre-truncation
	for i in *"_temp_"*".fastq"; do
		temp_fastq=`echo $temp_fastq","$i` 
	done
	temp_fastq="[""${temp_fastq:1}""]"

	python $hictoolPath"HiCtool_pre_truncation.py" -i $temp_fastq -e $restrictionEnzyme -p $threads
	
	cat "${fastq1%%.*}_temp_"*".trunc.fastq" > "${fastq1%%.*}.trunc.fastq"
	cat "${fastq2%%.*}_temp_"*".trunc.fastq" > "${fastq2%%.*}.trunc.fastq"

	rm *"_temp_"*".fastq"

	cat "${fastq1%%.*}"*"log"* > "${fastq1%%.*}_log.txt"
	cat "${fastq2%%.*}"*"log"* > "${fastq2%%.*}_log.txt"

	tot_reads_1=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_reads_2=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq2%%.*}_log.txt")

	tot_trunc_reads_1=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_trunc_reads_2=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq2%%.*}_log.txt")

	perc1=$(awk -v n1=$tot_trunc_reads_1 -v ntot1=$tot_reads_1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
	perc2=$(awk -v n1=$tot_trunc_reads_2 -v ntot1=$tot_reads_2 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)

	read_length=$(awk 'FNR == 1 {print $2}' "${fastq1%%.*}_log.txt")

	printf $fastq1"\n"$tot_reads_1" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_1" ("$perc1"%%)\n\n" > pre_truncation_log.txt
	printf $fastq2"\n"$tot_reads_2" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_2" ("$perc2"%%)" >> pre_truncation_log.txt

	rm *"_temp_"*
	rm "${fastq1%%.*}_log.txt"
	rm "${fastq2%%.*}_log.txt"
	
fi

fastq1_trunc="${fastq1%%.*}.trunc.fastq"
fastq2_trunc="${fastq2%%.*}.trunc.fastq"

echo -n "Aligning "$fastq1_trunc" ... "
(bowtie2 -p $threads -x $genomeIndex $fastq1_trunc -S HiCfile1.sam) 2>HiCfile1_log.txt
echo "Done!"
echo -n "Aligning "$fastq2_trunc" ... "
(bowtie2 -p $threads -x $genomeIndex $fastq2_trunc -S HiCfile2.sam) 2>HiCfile2_log.txt
echo "Done!"

# extracting the headers and read filtering
echo -n "Filtering HiCfile1.sam ... "
samtools view -H HiCfile1.sam > header1.txt
samtools view -F 4 -q 30 HiCfile1.sam > HiCfile1_hq.sam
echo "Done!"
echo -n "Filtering HiCfile2.sam ... "
samtools view -H HiCfile2.sam > header2.txt
samtools view -F 4 -q 30 HiCfile2.sam > HiCfile2_hq.sam
echo "Done!"

echo -n "Building log files ... "
n1=`wc -l HiCfile1_hq.sam | awk '{print $1}'`
n2=`wc -l HiCfile2_hq.sam | awk '{print $1}'`

nt1=`wc -l HiCfile1.sam | awk '{print $1}'`
h1=`wc -l header1.txt | awk '{print $1}'`
ntot1=`expr $nt1 - $h1`

nt2=`wc -l HiCfile2.sam | awk '{print $1}'`
h2=`wc -l header2.txt | awk '{print $1}'`
ntot2=`expr $nt2 - $h2`

perc1=$(awk -v n1=$n1 -v ntot1=$ntot1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
perc2=$(awk -v n2=$n2 -v ntot2=$ntot2 'BEGIN { print 100*n2/ntot2 }' | cut -c1-5)

printf "\n----------\n"$ntot1" reads; of these:\n  "$n1" ("$perc1"%%) aligned with MAPQ>=30" >> HiCfile1_log.txt
printf "\n----------\n"$ntot2" reads; of these:\n  "$n2" ("$perc2"%%) aligned with MAPQ>=30" >> HiCfile2_log.txt
echo "Done!"

rm HiCfile1.sam
rm HiCfile2.sam

echo -n "Selecting reads that are paired ... "
awk '{print $1}' HiCfile1_hq.sam | sort > readnames1.txt
awk '{print $1}' HiCfile2_hq.sam | sort > readnames2.txt
comm -12 readnames1.txt readnames2.txt > paired_reads.txt

if [ -z $max_lines ]
then
	# Select reads of the first sam file that are paires with the second sam file
	grep -Fwf paired_reads.txt HiCfile1_hq.sam | \
	cat header1.txt - | \
	samtools view -b -@ $threads - > HiCfile_pair1.bam
	rm HiCfile1_hq.sam

	# Select reads of the second sam file that are paired with the first sam file
	grep -Fwf paired_reads.txt HiCfile2_hq.sam | \
	cat header2.txt - | \
	samtools view -b -@ $threads - > HiCfile_pair2.bam
	rm HiCfile2_hq.sam

elif ! [ -z $max_lines ] && [ $max_lines -ge $fastq_lines ]
then
	# Select reads of the first sam file that are paires with the second sam file
	grep -Fwf paired_reads.txt HiCfile1_hq.sam | \
	cat header1.txt - | \
	samtools view -b -@ $threads - > HiCfile_pair1.bam
	rm HiCfile1_hq.sam

	# Select reads of the second sam file that are paired with the first sam file
	grep -Fwf paired_reads.txt HiCfile2_hq.sam | \
	cat header2.txt - | \
	samtools view -b -@ $threads - > HiCfile_pair2.bam
	rm HiCfile2_hq.sam

elif ! [ -z $max_lines ] && [ $max_lines -lt $fastq_lines ]
then

	# Splitting the paired reads file
	paired_reads=paired_reads.txt
	paired_reads_lines=`wc -l $paired_reads | awk '{print $1}'`
	max_lines_paired_reads=`expr $max_lines / 5`
	k=$max_lines_paired_reads
	count=1
	while [ $k -lt $paired_reads_lines ]
	do
		start=`expr $k - $max_lines_paired_reads + 1`
		quit=`expr $k + 1`
		sed -n "$start,"$k"p;"$quit"q" $paired_reads > "paired_reads_temp_"$count".txt"
		count=`expr $count + 1`
		k=`expr $k + $max_lines_paired_reads`
	done
	start=`expr $k - $max_lines_paired_reads + 1`
	sed -n "$start,"$paired_reads_lines"p" $paired_reads > "paired_reads_temp_"$count".txt"

	# Search for paired reads from each temporary file
	for i in "paired_reads_temp"*; do
		grep -Fwf $i HiCfile1_hq.sam > "HiCfile1_${i%%.*}.sam"
		grep -Fwf $i HiCfile2_hq.sam > "HiCfile2_${i%%.*}.sam"
	done

	cat "HiCfile1_paired_reads_temp"* > HiCfile1_paired.sam
	cat "HiCfile2_paired_reads_temp"* > HiCfile2_paired.sam

	rm *"temp"*

	cat header1.txt HiCfile1_paired.sam | \
	samtools view -b -@ $threads - > HiCfile_pair1.bam
	rm HiCfile1_paired.sam

	cat header2.txt HiCfile2_paired.sam | \
	samtools view -b -@ $threads - > HiCfile_pair2.bam
	rm HiCfile2_paired.sam

	rm HiCfile1_hq.sam
	rm HiCfile2_hq.sam

fi

echo "Done!"

echo -n "Updating log files ... "
n=`wc -l paired_reads.txt | awk '{print $1}'`

ntot1=`wc -l readnames1.txt | awk '{print $1}'`
ntot2=`wc -l readnames2.txt | awk '{print $1}'`

perc1=$(awk -v n1=$n -v ntot1=$ntot1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
perc2=$(awk -v n2=$n -v ntot2=$ntot2 'BEGIN { print 100*n2/ntot2 }' | cut -c1-5)

printf "; of these:\n    "$n" ("$perc1"%%) were paired and saved into HiCfile_pair1.bam" >> HiCfile1_log.txt
printf "; of these:\n    "$n" ("$perc2"%%) were paired and saved into HiCfile_pair2.bam" >> HiCfile2_log.txt
echo "Done!"

rm header1.txt
rm header2.txt
rm readnames1.txt
rm readnames2.txt
rm paired_reads.txt

echo "End data preprocessing: $(date)"
