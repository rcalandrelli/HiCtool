# Data preprocessing

This is the first section of the pipeline and it allows to pre-process the raw Hi-C data (fastq files), in order to generate input files for the normalization step. For more information about the Python functions used here check the [API documentation](https://sysbio.ucsd.edu/public/rcalandrelli/HiCtool_API_documentation.pdf).

## Table of Contents

1. [Preprocessing the data](#1-preprocessing-the-data)
   - [1.1. Downloading the raw data from GEO](#11-downloading-the-raw-data-from-geo)
2. [Creating the fragment-end (FEND) bed file](#2-creating-the-fragment-end-fend-bed-file)


## 1. Preprocessing the data

In order to start the Hi-C data analysis and preprocess your data you should have two fastq files, respectively of the first and the second reads of the pairs. If you wish to use public datasets on GEO and you need instructions to download the data, see [this section](#11-downloading-the-raw-data-from-geo). 

**Note!** To produce our final results, use this GEO accession number: **[GSM1551550](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1551550)**.

HiCtool allows to process Hi-C data generated using one or more of the following restriction enzymes:

- HindIII
- MboI
- DpnII
- Sau3AI
- BglII
- NcoI
- Hinfl

The **Arima Kit** uses a cocktail of restriction enzymes which includes MboI and Hinfl.
If your experiment was performed using a different restriction enzyme, please contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>.

In addition, the Bowtie2 genome index of the species under analysis should be provided too. If you do not have it, please run the following in order to generate it:
```unix
bowtie2-build hg38.fa index
```
```hg38.fa``` is the reference sequence in FASTA format (in this case for hg38), the output files in ``bt2`` format are named with the prefix ``index``.

The data preprocessing is performed with a single command line (replace parameters) and comprises the following steps:

- Pre-truncation of the reads that contain potential ligation junctions to keep the longest piece without a junction sequence ([Ay et al., 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7)).
- Independent mapping of the read pairs to the reference genome to avoid any proximity constraint.
- Removing the unmapped reads and selecting reads that were uniquely mapped with a MAPQ >= 30, i.e. the estimated probability of mapping error is <= 0.1%.

```unix
# Make the bash script executable
chmod u+x /HiCtool-master/scripts/HiCtool_run_preprocessing.sh

# Run the script
/HiCtool-master/scripts/HiCtool_run_preprocessing.sh \
-h /HiCtool-master/scripts/ \
-o your_output_directory \
-1 /myfastq_path/file1.fastq \
-2 /myfastq_path/file2.fastq \
-e MboI \
-g /path_to_the_genome_indexes/index \
-p 32
```
where:

- ``-h``: path where are the HiCtool scripts with the final trailing slash.
- ``-o``: path to save the output files. If the folder does not exist, it is created automatically.
- ``-1``: the fastq file with the first reads of the pairs.
- ``-2``: the fastq file with the second reads of the pairs.
- ``-e``: the restriction enzyme or enzymes passed between square brackets (example: [MboI,Hinfl] for the cocktail of the Arima Kit).
- ``-g``: Bowtie2 genome indexes. Only the filename should be passed here without extension.
- ``-p``: the number of parallel threads (processors) to use for alignment and preprocessing. The more the fastest the process.

The following output files are generated:

- ``HiCfile_pair1.bam`` and ``HiCfile_pair2.bam`` that are the bam files generated after alignment and filtering of the first and second reads in the pairs respectively.
- ``pre_truncation_log.txt`` with the information about the percentage of reads that have been truncated. This is also printed on the console:
```unix
SRR1658570_1.fastq
202095066 reads (length = 101 bp); of these:
29851195 (14.78%) contained a potential ligation junction and have been truncated.
SRR1658570_2.fastq
202095066 reads (length = 101 bp); of these:
28681691 (14.2%) contained a potential ligation junction and have been truncated.
```
The length distribution of the truncated reads is also plotted and saved to file.

![](/figures/SRR1658570_1.fastq_truncated_reads.png)

![](/figures/SRR1658570_2.fastq_truncated_reads.png)

- ``HiCfile1_log.txt`` and ``HiCfile2_log.txt`` are the log files with alignment and filtering statistics for the first and second reads in the pairs respectively.
```unix
HiCfile1_log.txt
202095066 reads; of these:
202095066 (100.00%) were unpaired; of these:
5770798 (2.86%) aligned 0 times
156759009 (77.57%) aligned exactly 1 time
39565259 (19.58%) aligned >1 times
97.14% overall alignment rate

----------
202095066 reads; of these:
172973813 (85.59%) aligned with MAPQ>=30; of these:
143415284 (82.91%) were paired and saved into HiCfile_pair1.bam

HiCfile2_log.txt
202095066 reads; of these:
202095066 (100.00%) were unpaired; of these:
13381441 (6.62%) aligned 0 times
149852422 (74.15%) aligned exactly 1 time
38861203 (19.23%) aligned >1 times
93.38% overall alignment rate

----------
202095066 reads; of these:
161438783 (79.88%) aligned with MAPQ>=30; of these:
143415284 (88.83%) were paired and saved into HiCfile_pair2.bam
```

### 1.1. Downloading the raw data from GEO

**Note!** If you have your fastq files generated from your custom experiment, you can skip this paragraph.

The source data in sra format are downloaded via GEO accession number using the command ``fastq-dump`` of [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump).

Before proceeding, you may need to [setup the output directory](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration) where the sra files will be saved. After having installed SRA Toolkit, go to the path where the software has been installed, under the subfolder “bin”, and run the following command line:
```unix
./vdb-config -i
```
This will open an interface that will allow you to setup/change your output directory.

To download the data related to a GEO accession number, go to the bottom of that page and click on the SRA number under the section “Relations”. After that, under the section “Runs” you will find the SRR files, then run the following:
```unix
fastq-dump SRRXXXXXXX --split-3
```
where ``SRRXXXXXXX`` has to be replaced with the specific number of the run you want to download (**[SRR1658570](https://www.ncbi.nlm.nih.gov/sra?term=SRX764936) in this documentation**).

To be more specific, this code will either download the SRA file under the output directory you have set with the interface above, but it will also convert the SRA file into fastq files and dump each read into separate files in the current working directory, to produce:

- SRRXXXXXXX_1.fastq
- SRRXXXXXXX_2.fastq
- SRRXXXXXXX.fastq

where paired-end reads in SRRXXXXXXX.sra are split and stored into **SRRXXXXXXX_1.fastq** and **SRRXXXXXXX_2.fastq**, and SRRXXXXXXX.fastq (if present) contains reads with no mates.

**Note!** To produce our final results, use this GEO accession number: **[GSM1551550](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1551550)**.


## 2. Creating the fragment-end (FEND) bed file

The fragment-end (FEND) bed file is used to normalize the data and it contains restriction site coordinates and additional information related to fragment properties (GC content and mappability score). Specifically, for each fragment the GC content of 200 bp upstream and downstream to the restriction site is computed. For the mappability score, the entire genome sequence is split into artificial reads (50 bp reads, starting every 10 bp) and then mapped back to the genome. For each fragment end the mappability score is then defined to be the portion of artificial reads mapped uniquely to the genome (MAPQ > 30) within a 500-bp window upstream and downstream to the fragment. Fragment ends with a mappability score less than 0.5 are then discarded (Yaffe and Tanay, 2011).

Since the following steps may be time consuming, we provide the most common FEND files available for download (DpnII is the same restriction site than MboI):

- [HindIII-hg38](http://data.genomegitar.org/HindIII_hg38_gc_map_valid.zip)
- [MboI-hg38](http://data.genomegitar.org/MboI_hg38_gc_map_valid.zip) (file used in this documentation)
- [NcoI-hg38](http://data.genomegitar.org/NcoI_hg38_gc_map_valid.zip)
- [HindIII-mm10](http://data.genomegitar.org/HindIII_mm10_gc_map_valid.zip)
- [MboI-mm10](http://data.genomegitar.org/MboI_mm10_gc_map_valid.zip)

**Perform the following steps ONLY if you need to generate a new fragment end bed file (because you are using another species or a different restriction enzyme than those provided above). Otherwise, download the file of interest and go to the [data normalization section](https://github.com/Zhong-Lab-UCSD/HiCtool/tree/master/tutorial#2-data-normalization-and-visualization).**

***

In order to align all the restriction sites for a certain cutting enzyme, a ``fastq`` file related to the enzyme cutting site has to be provided. For the quality score of the restriction enzyme sequence, we can simply add a default average score ``I``:
```unix
echo -e "@HindIII\nAAGCTT\n+\nIIIIII" > HindIII.fastq
echo -e "@MboI\nGATC\n+\nIIII" > MboI.fastq
echo -e "@NcoI\nCCATGG\n+\nIIIIII" > NcoI.fastq
```
After this, implement the multiple alignment command in Bowtie 2 to locate all the coordinates of the restriction enzyme sites:
```unix
(bowtie2 -p 32 -k 10000000 -x your_genome_index -U MboI.fastq -S restrictionsites.sam) 2>restrictionsites_log.txt
```
where the ``-k`` argument changes Bowtie 2 research behavior. By default, ``bowtie2`` searches for distinct, valid alignments for each read. When it finds a valid alignment, it continues looking for alignments that are nearly as good or better and the best alignment found is reported. When ``-k <int>`` is specified, ``bowtie2`` searches for at most ``<int>`` distinct, valid alignments for each read. The search terminates when it can not find more distinct valid alignments, or when it finds ``<int>``, which happens first.

In order to use the restriction sites file as an input for the following analysis, we have to convert the ``sam`` file to ``bed`` file via [SAMtools](http://samtools.sourceforge.net/) and [bedtools](http://bedtools.readthedocs.org/en/latest/):
```unix
samtools view -b restrictionsites.sam | bedtools bamtobed -i > restrictionsites.bed
```
- **If you will be using the [Hi-Corrector normalization approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-matrix-balancing.md), you do NOT need to perform the following steps, because additional information such as GC content or mappability score are not needed. You can use ``restrictionsites.bed`` as the FEND file.**
- **Proceed with the following steps if you will be using [Yaffe and Tanay's normalization approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-yaffe-tanay.md).**

***

First ``restrictionsites.bed`` is split into separate files, one per each chromosome (**update the chromosomes list according to your species**):
```unix
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for i in "${chromosomes[@]}"; do
awk -v var="$i" '(NR>1) && ($1==var)' restrictionsites.bed > $i.bed
done
```

The GC content and mappability score information must be in comma-separated format. First, each feature is computed for each separate chromosome, finally the files are merged together and parsed to generate a unique bed file. For this part, **Python multiprocessing** is used to consistently reduce the computation time. It is recommended to use the highest number of threads available in your processor (the highest up to the total number of chromosomes of your species).

- Download the GC content information for the species of interest from the UCSC website at http://hgdownload.cse.ucsc.edu/gbdb/hg38/bbi/ (replace the species name in the link if needed) and use the file ``gc5Base.bw``. This file is in BigWig format, before running step 2) we need to convert it to BedGraph format (tab separated file with 4 columns: chromosome, start, end, score; in this case "score" is the GC content). After this, the file has to be splitted into separate txt files, one per each chromosome, named as **chr1.txt, chr2.txt, … , chrX.txt, chrY.txt**. To do so, run the following Unix script (note that given the dimension of the files, this process may require sometime):
```unix
wget link_to_gc5Base.bw
bigWigToBedGraph gc5Base.bw gc5Base.bedGraph

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
for i in "${chromosomes[@]}"; do
awk -v var="$i" '(NR>1) && ($1==var)' gc5Base.bedGraph | awk -v OFS='\t' '{print $1, $2, $3, $4}' > $i.txt
done
```
- Add the GC content information using the Python script [HiCtool_add_fend_gc_content.py](/scripts/HiCtool_add_fend_gc_content.py). Open the script, **update the parameters on the top and save**. Then just execute the script to add the gc content information (using 24 threads, we took around 9 hours for all the chromosomes of hg38-MboI):
```Python
execfile('HiCtool_add_fend_gc_content.py')
```
- Generate artificial reads using the Python function inside [HiCtool_artificial_reads.py](/scripts/HiCtool_artificial_reads.py) (this step is required only once per reference genome):
```Python
execfile('HiCtool_artificial_reads.py')
generate_artificial_reads(genome_file, output_reads_file)
```
where ``genome_file`` is the reference genome sequence in ``fasta`` format, ``output_reads_file`` is the file where to save reads in ``fastq`` format.
- Align the artificial reads to the reference genome, remove unmapped reads and generate a bed formatted file where the last field is the MAPQ score (this step is required only once per reference genome). This is the same format of the GC content files downloaded from UCSC, where the last field was the GC percentage instead of MAPQ:
```unix
(bowtie2 -p 32 -x your_genome_index artificial_reads.fastq -S artificial_reads.sam) 2>artificial_reads.log
samtools view -F 4 artificial_reads.sam > artificial_reads_mapped.sam
awk -v OFS='\t' '{print $3, $4-1, $4-1+50, $5}' artificial_reads_mapped.sam > artificial_reads_mapped.txt
```
- Split the mapped reads into separate files, one per chromosome (**update the list of chromosome names if a different species is used**):
```unix
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
for i in "${chromosomes[@]}"; do
awk -v var="$i" '(NR>1) && ($1==var)' artificial_reads_mapped.txt | awk -v OFS='\t' '{print $1, $2, $3, $4}' > $i.txt
done
```
- Add the mappability score information using the Python script [HiCtool_add_fend_mappability.py](/scripts/HiCtool_add_fend_mappability.py). Open the script, **update the parameters on the top and save**. Then just execute the script to add the mappability information (using 24 threads, we took around 9 hours for all the chromosomes of hg38-MboI):
```Python
execfile('HiCtool_add_fend_mappability.py')
```
- Sort by coordinate, merge the files together, remove fragment ends with a mappability score < 0.5, parse GC content and mappability score in comma separated format and add the header:
```unix
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for i in "${chromosomes[@]}"; do
sort -k 2,2n "$i"_restrictionsites_gc_map.bed | cat >> restrictionsites_gc_map.bed
done

awk '(NR>1) && ($9 >= 0.5) && ($10 >= 0.5)' restrictionsites_gc_map.bed | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7 "," $8, $9 "," $10}' > restrictionsites_gc_map_valid_noHeader.bed
echo -e 'chr\tstart\tstop\tname\tscore\tstrand\tgc\tmappability' | cat - restrictionsites_gc_map_valid_noHeader.bed > restrictionsites_gc_map_valid.bed

rm restrictionsites_gc_map_valid_noHeader.bed
```
**restrictionsites_gc_map_valid.bed** is the final FEND bed file that will be used in the normalization pipeline to remove biases.
