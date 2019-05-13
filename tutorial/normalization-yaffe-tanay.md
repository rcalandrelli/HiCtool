# Data normalization with explicit-factor correction model of Yaffe and Tanay

This pipeline illustrates the procedure to normalize and visualize Hi-C **intra-chromosomal contact data only** following the explicit-factor model of [Yaffe and Tanay](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html).

For more information about the Python functions used here check the [API documentation](https://sysbio.ucsd.edu/public/rcalandrelli/HiCtool_API_documentation.pdf).

## Table of contents

1. [Running HiFive functions](#1-running-hifive-functions)
2. [Normalizing the data](#2-normalizing-the-data)
   - [2.1. Normalized fend data](#21-normalized-fend-data)
   - [2.2. Normalized enrichment data](#22-normalized-enrichment-data)
   - [2.3. Multi-processing normalization](#23-multi-processing-normalization)
3. [Visualizing the data](#3-visualizing-the-data)
   - [3.1. Visualizing the contact data](#31-visualizing-the-contact-data)
   - [3.2. Visualizing the enrichment data](#32-visualizing-the-enrichment-data)

## 1. Running HiFive functions

We resort to the HiFive package, and specifically the binning algorithm, to normalize the data using the approach of Yaffe and Tanay. HiFive allows to remove spurious ligation products, as well as PCR duplicates and non-informative reads before generating the contact matrix.

The Python script [HiCtool_hifive.py](/scripts/HiCtool_hifive.py) is used to run all the steps needed in order to generate the data used for normalization: 

- Creating the Fend object.
- Creating the HiCData object from a Fend object and mapped data in bam format. At this step spurious ligation products (paired-reads whose total distance to their respective restriction sites exceeds 500 bp) are removed. In addition, PCR duplicates are removed and reads with ends mapping to the same fragment and reads with ends mapping to adjacent fragments on opposite strands are also excluded, to consider the possibility of incomplete restriction enzyme digestion and fragment circularization.
- Creating the HiC project object, which stores the information about the observed contact data that we will use for the downstream analysis.
- Filtering fragments that do not have at least one interaction before learning correction parameters.
- Estimating the distance-dependence relationship from the data prior to normalization, in order to avoid biases that may result due to restriction site distribution characteristics or the influence of distance/signal relationship. Restriction sites over the genome are unevenly distributed and this results in a large set of distances between fragments and their neighbors. Since the interaction frequency is strongly inversely-related to inter-fragment distance, this means that fragments surrounded by shorter ones will show higher nearby interactions than those with longer adjacent fragments, due to the uneven distribution of the restriction sites position.
- Learning the correction model for Hi-C data. For the normalization, we take into account of fragments length, inter-fragment distance, GC content and mappability score biases, according to the information included in the Fend object. We also consider a minimum distance of 500 kb between fragments to take into account of the effect of biological biases (TSSs and CTCF bound sites) while learning the correction parameters.

For more information about these functions, please see [HiFive’s API documentation](http://bxlab-hifive.readthedocs.org/en/latest/api.html). To run these steps execute the following command on the Unix console (update parameters properly):
```unix
python /HiCtool-master/scripts/HiCtool_hifive.py \
-f restrictionsites.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e MboI \
-m Yaffe-Tanay
```
where:

- ``-f`` is the FEND file in bed format from preprocessing.
- ``--b1`` is the first bam file from preprocessing.
- ``--b2`` is the second bam file from preprocessing.
- ``-e`` is the restriction enzyme or enzymes names between square brackets (example [MboI,Hinfl]).
- ``-m`` is the normalization model used (Yaffe_Tanay in this case).

**The following output files are generated:**

- ``fend_object.hdf5``
- ``HiC_data_object.hdf5``
- ``HiC_project_object.hdf5``
- ``HiC_norm_binning.hdf5`` to be used in the following section.


## 2. Normalizing the data

For the normalization, observed data and correction parameters to remove biases to obtain the corrected read counts are required. Therefore, the observed contact matrix and the fend expected contact matrix are calculated. In addition, the enrichment expected contact matrix is calculated to compute the observed over expected enrichment values, considering also the distance between fends.

For each chromosome, the following five matrices are computed at a bin size of 40 kb (the bin size can be changed with a function parameter). Every contact matrix is AUTOMATICALLY saved in txt format using the function ``save_matrix``.

Data are compressed in a format to reduce storage occupation and improving saving and loading time ([see here for more details](/tutorial/HiCtool_compressed_format.md)).

- The **observed data** contain the observed reads count for each bin.
- The **fend expected data** contain the learned correction value to remove biases related to fends for each bin.
- The **enrichment expected data** contain the expected reads count for each bin, considering the linear distance between read pairs and the learned correction parameters.
- The **normalized fend data** contain the corrected read count for each bin.
- The **normalized enrichment data** ("observed over expected" matrix) contain the enrichment value (O/E) for each bin.

First execute the script [HiCtool_yaffe_tanay.py](/scripts/HiCtool_yaffe_tanay.py) in the Python or iPython console:
```Python
execfile('HiCtool_yaffe_tanay.py')
```
### 2.1. Normalized fend data

To calculate and save the **normalized intra-chromosomal contact matrix** for a chromosome ``a_chr``, use the function ``normalize_chromosome_fend_data``:
```Python
fend_normalized_chr6 = normalize_chromosome_fend_data(a_chr='6', bin_size=40000, 
                                                      input_file='HiC_norm_binning.hdf5', 
                                                      species='hg38',
                                                      save_obs=True, save_expect=False)
```
Data are compressed in a format to reduce storage occupation and improving saving and loading time ([see here for more details](/tutorial/HiCtool_compressed_format.md)). To load a previously generated contact matrix use the function ```load_matrix```:
```Python
my_contact_matrix = load_matrix('my_contact_matrix.txt')
```
where ``'my_contact_matrix.txt'`` is a contact matrix file saved using ``normalize_chromosome_fend_data`` .

### 2.2. Normalized enrichment data

To calculate and save the **"observed/expected" intra-chromosomal contact matrix** for a chromosome ``a_chr`` use the function ``normalize_chromosome_enrich_data`` (see API Documentation):
```Python
enrich_normalized_chr6 = normalize_chromosome_enrich_data(a_chr='6', bin_size=40000, 
                                                          input_file='HiC_norm_binning.hdf5', 
                                                          species='hg38',
                                                          save_obs=True, save_expect=False)
```
**Note!**
If you need only the normalized contact matrices, there is no need to calculate also the enrichment data. If you do not need the expected data, do not save it since they are the biggest files and the process may take time.

Data are compressed in a format to reduce storage occupation and improving saving and loading time ([see here for more details](/tutorial/HiCtool_compressed_format.md)). To load a previously generated contact matrix use the function ```load_matrix```:
```Python
my_contact_matrix = load_matrix('my_contact_matrix.txt')
```
where ``'my_contact_matrix.txt'`` is a contact matrix file saved using ``normalize_chromosome_enrich_data``.

### 2.3. Multi-processing normalization

To calculate and save the normalized contact matrices in parallel for multiple chromosomes, use the script [HiCtool_normalize_fend_parallel.py](/scripts/HiCtool_normalize_fend_parallel.py). **Open the script, update the parameters on the top and save.** Then, just execute the script:
```Python
execfile('HiCtool_normalize_fend_parallel.py')
```
To calculate and save the "observed/expected" contact matrices in parallel use the script [HiCtool_normalize_enrich_parallel.py](/scripts/HiCtool_normalize_enrich_parallel.py). **Open the script, update the parameters on the top and save.** Then, just execute the script:
```Python
execfile('HiCtool_normalize_enrich_parallel.py')
```
## 3. Visualizing the data

To plot the contact maps, first execute the script [HiCtool_yaffe_tanay.py](/tutorial/HiCtool_yaffe_tanay.py) in the Python or iPython console:
```Python
execfile('HiCtool_yaffe_tanay.py')
```
### 3.1. Visualizing the contact data

This part is to plot heatmaps and histograms of the contact data. 

To plot and save the heatmap and histogram use the function ```plot_chromosome_data```:
```Python
plot_chromosome_data('HiCtool_chr6_40kb_normalized_fend.txt', 
                     a_chr='6', bin_size=40000, full_matrix=False, 
                     start_coord=50000000, end_coord=54000000, 
                     species='hg38', 
                     data_type="normalized_fend", 
                     my_colormap=['white', 'red'], 
                     cutoff_type='percentile', cutoff=95, max_color='#460000', 
                     my_dpi=1000, 
                     plot_histogram=True)
```
Instead of ``'HiCtool_chr6_40kb_normalized_fend.txt'``, the object containing the contact matrix calculated above ``fend_normalized_chr6`` can be passed as well.

This function can be used also to plot observed data, expected fend and enrichment data by simply passing a different input matrix as first parameter.

Heatmap             |  Histogram
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_40kb_normalized_fend.png)  |  ![](/figures/HiCtool_chr6_40kb_normalized_fend_histogram.png)


**Additional example of the contact matrix for chromosome 6 at 1 Mb resolution**

In order to change the heatmap resolution, first data have to be normalized at the desired resolution set with the parameter ``bin_size`` of ``normalize_chromosome_fend_data`` ([see section 2.1.](#21-normalized-fend-data)):
```Python
fend_normalized_chr6 = normalize_chromosome_fend_data(a_chr='6', 
                                                      bin_size=1000000, 
                                                      input_file='HiC_norm_binning.hdf5', 
                                                      species='hg38',
                                                      save_obs=True, 
                                                      save_expect=False)
```
Then, we plot the entire heatmap (we also change here the color map to white and blue):
```Python
plot_chromosome_data(fend_normalized_chr6, 
                     a_chr='6', 
                     bin_size=1000000, 
                     full_matrix=True, 
                     species='hg38', 
                     data_type="normalized_fend", 
                     my_colormap=['white', 'blue'], 
                     cutoff_type='percentile', cutoff=95, max_color='#460000', 
                     my_dpi=1000, 
                     plot_histogram=False)
```
![](/figures/HiCtool_chr6_1mb_normalized_fend.png)

### 3.2. Visualizing the enrichment data

This part is to plot the heatmap and histogram for the enrichment normalized data ("observed over expected"). The **log2 of the data** is plotted to quantify the positive enrichment (red) and the negative enrichment (blue). Loci (pixels) equal to zero before performing the log2 (deriving from zero observed contacts) are shown in gray. Loci (pixels) where enrichment expected contact was zero before performing the ratio (observed / expected) are shown in black.

To plot and save the heatmap and histogram use the function ```plot_chromosome_enrich_data```:
```Python
plot_chromosome_enrich_data('HiCtool_chr6_40kb_normalized_enrich.txt', 
                            a_chr='6', 
                            bin_size=40000, 
                            full_matrix=False, 
                            start_coord=50000000, end_coord=54000000, 
                            species='hg38', 
                            my_dpi=1000, 
                            plot_histogram=True)
```
Instead of ``'HiCtool_chr6_40kb_normalized_enrich.txt'``, the object containing the contact matrix calculated above ``enrich_normalized_chr6`` can be passed as well.

Heatmap             |  Histogram
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_40kb_normalized_enrich.png)  |  ![](/figures/HiCtool_chr6_40kb_normalized_enrich_histogram.png)

**Additional example of the enrichment contact matrix for chromosome 6 at 1 Mb resolution**

In order to change the heatmap resolution, first data have to be calculated at the desired resolution set with the parameter ``bin_size`` of ``normalize_chromosome_enrich_data`` ([see section 2.2.](#22-normalized-enrichment-data)):
```Python
enrich_normalized_chr6 = normalize_chromosome_enrich_data(a_chr='6', 
                                                          bin_size=1000000, 
                                                          input_file='HiC_norm_binning.hdf5', 
                                                          species='hg38',
                                                          save_obs=False, 
                                                          save_expect=False)
```
Then, we plot the entire heatmap with a maximum and minimum cutoff for the log2 at 4 and -4 respectively:
```Python
plot_chromosome_enrich_data(enrich_normalized_chr6, 
                            a_chr='6', 
                            bin_size=1000000, 
                            full_matrix=True, 
                            species='hg38',
                            cutoff_max=4,
                            cutoff_min=-4,
                            plot_histogram=True)
```
Heatmap             |  Histogram
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_1mb_normalized_enrich.png)  |  ![](/figures/HiCtool_chr6_1mb_normalized_enrich_histogram.png)
