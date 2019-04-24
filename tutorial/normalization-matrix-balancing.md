# Data normalization with the matrix balancing approach of Hi-Corrector

This pipeline illustrates the procedure to normalize and normalize a **global Hi-C contact map** (intra- and inter-chromosomal interactions) following the matrix balancing approach of [Hi-Corrector](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html).

For more information about the Python functions used here check the [API documentation](https://sysbio.ucsd.edu/public/rcalandrelli/HiCtool_API_documentation.pdf).

## Table of Contents

1. [Running HiFive functions](#1-running-hifive-functions)
2. [Generating the global observed contact matrix](#2-generating-the-global-observed-contact-matrix)
3. [Normalizing the global contact matrix](#3-normalizing-the-global-contact-matrix)
4. [Visualizing the data](#4-visualizing-the-data)
   - [4.1. Visualizing the global contact data](#41-visualizing-the-global-contact-data)
   - [4.2. Visualizing a single heatmap](#42-visualizing-a-single-heatmap)

## 1. Running HiFive functions

We resort to the HiFive package in order to generate the global observed contact matrix (intra- and inter-chromosomal contact maps all together to form a single, global contact matrix). HiFive allows to remove spurious ligation products, as well as PCR duplicates and non-informative reads.

The Python script [HiCtool_hifive.py](/scripts/HiCtool_hifive.py) is used to run all the three steps needed in order to obtain the observed contact data, whose outputs are ``.hdf5`` files: 

- Creating the Fend object.
- Creating the HiCData object from a Fend object and mapped data in bam format. At this step spurious ligation products (paired-reads whose total distance to their respective restriction sites exceeds 500 bp) are removed. In addition, PCR duplicates are removed and reads with ends mapping to the same fragment and reads with ends mapping to adjacent fragments on opposite strands are also excluded, to consider the possibility of incomplete restriction enzyme digestion and fragment circularization.
- Creating the HiC project object, which stores the information about the observed contact data that we will use for the downstream analysis.

For more information about these functions, please see [HiFive’s API documentation](http://bxlab-hifive.readthedocs.org/en/latest/api.html). To run these steps execute the following command on the Unix console (update parameters properly):
```unix
python /HiCtool-master/scripts/HiCtool_hifive.py \
-f restrictionsites.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e MboI \
-m Hi-Corrector
```
where:

- ``-f`` is the FEND file from preprocessing.
- ``--b1`` is the first bam file from preprocessing.
- ``--b2`` is the second bam file from preprocessing.
- ``-e`` is the restriction enzyme or enzymes names between square brackets.
- ``-m`` is the normalization model used (Yaffe-Tanay or Hi-Corrector).

**The following output files are generated:**

- ``fend_object.hdf5``
- ``HiC_data_object.hdf5``
- ``HiC_project_object.hdf5`` to be used in the following section.


## 2. Generating the global observed contact matrix

This section will allow you to generate a **global square observed contact matrix** (24-by-24 chromosomes for hg38). The total number of bins of this big matrix will depend on the resolution of the data, and it can be estimated as the entire genome length over the resolution (for hg38 at 1Mb resolution is 3078x3078). The observed data contain the observed read count per each bin.

Especially at higher resolution, the generation of the global observed contact matrix may be computationally expensive and require long time. Therefore, we implemented a code to allow job parallelization. Each row of the contact matrix is computed in parallel, meaning all the contact matrices per each chromosome, and finally they are merged together to generate the global matrix. Each row of the matrix is saved in a temporary file, which is automatically deleted after the job is done.

To calculate and save the global observed contact matrix use the script [HiCtool_global_map_observed.py](/scripts/HiCtool_global_map_observed.py) and run this command:
```unix
python /HiCtool-master/scripts/HiCtool_global_map_observed.py \
-i HiCtool_project_object.hdf5 \
-o /output_path/ \
-b 1000000 \
-s hg38
-c /HiCtool-master/scripts/chromSizes/ \
--save_each 0 \
-p 24
```
where:

- ``-i``: Project object file in ``.hdf5`` format obtained with ``HiCtool_hifive.py``.
- ``-o``: Output path to save the observed contact matrix with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``--save_each``: Set to 1 to save each single contact matrix, 0 otherwise.
- ``-p``: Number of parallel threads to use. It has to be less or equal than the number of chromosomes of your species.

**The following output files are generated:**

- ``HiCtool_1mb_matrix_global_observed.txt``, the global matrix saved using a compressed format ([see here for more details](/tutorial/HiCtool_compressed_format.md)).
-``HiCtool_1mb_matrix_global_observed_tab.txt``, the global matrix saved in tab separated format. This matrix will be used in the next section to normalize the data.
- ``info.txt``, which contains the number of rows of the global matrix and the average rowsum.


After having generated the global observed contact matrix, it is possible to extract a single contact matrix (either intra- or inter-chromosomal) using the function ``extract_single_map`` of [HiCtool_full_map_analysis.py](/scripts/HiCtool_full_map_analysis.py) as following:
```Python
execfile('HiCtool_full_map_analysis.py')
global_observed = load_matrix('HiCtool_1mb_matrix_global_observed.txt')

chr1_intra = extract_single_map(input_global_matrix=global_observed, tab_sep=False, 
                                chr_row='1', chr_col='1', bin_size=1000000, 
                                data_type='observed', save_output=True, save_tab=True)

chr1_2_inter = extract_single_map(input_global_matrix=global_observed, tab_sep=False, 
                                  chr_row='1', chr_col='2', bin_size=1000000, 
                                  data_type='observed', save_output=True, save_tab=True)
```
**Tip!** ``extract_single_map`` can accept also directly the path to the global matrix file (``input_global_matrix='HiCtool_1mb_matrix_global_observed.txt'``) however, especially at higher resolution, the loading step of the matrix may require long time. Therefore, it is suggested to load once the matrix in the workspace using ``load_matrix`` and then work with it.


## 3. Normalizing the global contact matrix

Here we normalize the data using the sequential implementation from Hi-Corrector (ic_mes) which is "memory efficient and can run on any single computer with limited memory, even for Hi-C datasets of large size. It is designed to overcome the memory limit by loading a portion of data into the memory at each time." ([Li, Wenyuan, et al. "Hi-Corrector: a fast, scalable and memory-efficient package for normalizing large-scale Hi-C data." Bioinformatics 31.6 (2014): 960-962](https://academic.oup.com/bioinformatics/article/31/6/960/215261)).

The Hi-Corrector source code ([see here](https://github.com/Zhong-Lab-UCSD/HiCtool#installation)) is already inside ``/HiCtool-master/scripts/``. To normalize the data, run the following command:
```unix
# Make the bash script executable
chmod u+x /HiCtool-master/scripts/HiCtool_run_ic_mes.sh

# Run the script
/HiCtool-master/scripts/HiCtool_run_ic_mes.sh \
-q 100 \
-m 100 \
-r 3078 \
-s 17237 \
-h /HiCtool-master/scripts/Hi-Corrector1.2/ \
-i /output_path/HiCtool_1mb_matrix_global_observed_tab.txt
```
where:

- ``-q``: maximum number of iterations performed in the algorithm.
- ``-m``: the memory size. Its unit is Megabytes (MB).
- ``-r``: the number of rows or columns of the input chromatin contact frequency matrix to be normalized (provided in  ``info.txt`` generated in [section 2](#2-generating-the-global-observed-contact-matrix)).
- ``-s``: the row sum after normalization. The iterative correction algorithm can allow users to specify the row sum after the normalization, because this method is a matrix scaling approach that normalizes the matrix to be a doubly stochastic matrix (rows and columns sums equal to 1). Then we can multiple each element of this normalized matrix by the given value of this parameter, say 10.0 or 100.0 or whatever you choose. In such a way, the row sums of normalized matrix becomes this number (10.0 or 100.0 or whatever you choose). In ``info.txt`` we provide a row sum value that you could use calculated as "the average number of contacts of the observed matrix multiplied by the number of rows" to make the normalized data counts "comparable" with the observed ones. The choice is arbitrary.
- ``-h``: the path to the Hi-Corrector source code with the final trailing slash ``/``.
- ``-i``: the observed global contact matrix in tab delimited format.

This command creates a **folder named ``output_ic_mes``** with 3 files inside:

- ``output.log``: a log file
- ``output.bias``: a bias file used by the software to normalize the data
- ``output_normalized.txt``: the **global normalized contact matrix** in tab separated format

After having normalized the data, it is possible to extract a single normalized contact matrix (either intra- or inter-chromosomal) using the function ``extract_single_map`` of [HiCtool_full_map_analysis.py](/scripts/HiCtool_full_map_analysis.py) as following:
```Python
execfile('HiCtool_full_map_analysis.py')
global_normalized = load_matrix_tab("output_ic_mes/output_normalized.txt")

chr1_intra_norm = extract_single_map(input_global_matrix=global_normalized, tab_sep=True, 
                                     chr_row='1', chr_col='1', bin_size=1000000,       
                                     data_type='normalized', save_output=True, save_tab=True)  

chr1_2_inter_norm = extract_single_map(input_global_matrix=global_normalized, tab_sep=True, 
                                       chr_row='1', chr_col='2', bin_size=1000000, 
                                       data_type='normalized', save_output=True, save_tab=True)
```
**Tip!** ``extract_single_map`` can accept also directly the path to the global matrix file (``input_global_matrix='output_ic_mes/output_normalized.txt'``) however, especially at higher resolution, the loading step of the matrix may require long time. Therefore, it is suggested to load once the matrix in the workspace using ``load_matrix_tab`` and then work with it.

## 4. Visualizing the data

To plot the contact maps use the function ``plot_map`` of [HiCtool_full_map_analysis.py](/scripts/HiCtool_full_map_analysis.py).
```Python
execfile('HiCtool_full_map_analysis.py')
global_observed = load_matrix('HiCtool_1mb_matrix_global_observed.txt')
global_normalized = load_matrix_tab('output_ic_mes/output_normalized.txt')
```
**Tip!** ``plot_map`` used below can accept also directly the path to the global matrix files loaded above however, especially at higher resolution, the loading step of the matrix may require long time. Therefore, it is suggested to load once the matrix of interest in the workspace using ``load_matrix`` or ``load_matrix_tab`` as appropriate, and then plot.

### 4.1. Visualizing the global contact data

You can visualize either the observed or the normalized data. Here we plot both the global maps at 1 Mb resolution as calculated above.
```Python
# Observed data
plot_map(input_matrix=global_observed, isGlobal=True,
         bin_size=1000000, data_type='observed', species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc', cutoff=99, max_color='#460000')
```
![](/figures/HiCtool_1mb_observed.png)

```Python
# Normalized data
plot_map(input_matrix=global_normalized, isGlobal=True,
         bin_size=1000000, data_type='normalized', species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc', cutoff=99, max_color='#460000')
```
![](/figures/HiCtool_1mb_normalized.png)

### 4.2. Visualizing a single heatmap

A single contact matrix can be plotted by passing as argument the chromosome in the rows (``chr_row``) and in the columns (``chr_col``). 

To plot the **intra-chromosomal heatmap** of chromosome 6, run the following:
```Python
# Observed contact heatmap
plot_map(input_matrix=global_observed, isGlobal=True,
         chr_row='6', chr_col='6', bin_size=1000000, 
         data_type="observed", species='hg38',
         my_colormap=['white', 'red'], cutoff_type='perc',
         cutoff=99, max_color='#460000')

# Normalized contact heatmap
plot_map(input_matrix=global_normalized, isGlobal=True,
         chr_row='6', chr_col='6', bin_size=1000000, 
         data_type="normalized", species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc', cutoff=99, max_color='#460000')
```
Observed (chr 6)           |  Normalized (chr 6)
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_chr6_1mb_observed.png)  |  ![](/figures/HiCtool_chr6_chr6_1mb_normalized.png)

An **inter-chromosomal heatmap** can be also plotted (chr6-chr3) by setting the parameters ``chr_row`` and ``chr_col`` (we plot also the histogram of the contact distribution):
```Python
# Observed contact heatmap
plot_map(input_matrix=global_observed, isGlobal=True,
         chr_row='6', chr_col='3', bin_size=1000000, 
         data_type="observed", species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc', cutoff=99, max_color='#460000',
         plot_histogram=True)

# Normalized contact heatmap
plot_map(input_matrix=global_normalized, isGlobal=True,
         chr_row='6', chr_col='3', bin_size=1000000, 
         data_type="normalized", species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc', cutoff=99, max_color='#460000',
         plot_histogram=True)
```
Observed (chr6-chr3)            |  Normalized (chr6-chr3)
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_chr3_1mb_observed.png)  |  ![](/figures/HiCtool_chr6_chr3_1mb_normalized.png)
![](/figures/HiCtool_chr6_chr3_1mb_observed_histogram.png)  |  ![](/figures/HiCtool_chr6_chr3_1mb_normalized_histogram.png)

In addition, only a **region of the heatmap** can be plotted by setting the parameters ``chr_row_coord`` and ``chr_col_coord``. These are lists with two integers indicating the start and end coordinate of the chromosome on the rows and on the columns respectively.
```Python
# Intra-chromosomal map
plot_map(input_matrix=global_normalized, isGlobal=True,
         chr_row='6', chr_col='6', bin_size=1000000, 
         chr_row_coord=[0,80000000], chr_col_coord=[0,80000000],
         data_type="normalized", species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc', cutoff=99, max_color='#460000',
         plot_histogram=False)

# Inter-chromosomal map
plot_map(input_matrix=global_normalized, isGlobal=True,
         chr_row='6', chr_col='3', bin_size=1000000, 
         chr_row_coord=[0,50000000], chr_col_coord=[0,80000000],
         data_type="normalized", species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc', cutoff=99, max_color='#460000',
         plot_histogram=False)
```
Normalized (chr6) 0-80 Mb         |  Normalized (chr6-chr3) 0-50Mb; 0-80Mb
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_chr6_1mb_0-80mb_normalized.png)  |  ![](/figures/HiCtool_chr6_chr3_1mb_0-50_0-80mb_normalized.png)
