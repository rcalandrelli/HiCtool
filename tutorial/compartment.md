# A/B compartment analysis

This pipeline illustrates the procedure to calculate principal components (PC) of the Pearson correlation matrix that can be used to delineate A/B compartments in Hi-C data at low resolution (usually 1 mb or 500 kb). The code allows to calculate both PC1 and PC2. Usually, the sign of the eigenvector (PC1) indicates the compartment.

## Table of contents

1. [Calculating the principal component](#1-calculating-the-principal-component)
2. [Plotting the principal component](#2-plotting-the-principal-component)

## 1. Calculating the principal component

HiCtool allows to calculate either the first (typically used) or the second principal component of the Person correlation matrix. In order to do so, the Person correlation matrix has to be calculated first as presented [here](/tutorial/normalization-yaffe-tanay.md#22-normalizing-enrichment-data-and-calculating-the-person-correlation-matrix).

python /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_compartment_analysis.py \
--action calculate_pc \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1

where:

- ``--action``: action to perform (here ``normalize_fend``).
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--chr``: The chromosome to normalize.
- ``--save_obs``: Set to 1 to save the observed contact matrix, 0 otherwise.
- ``--save_expect``: Set to 1 to save the fend expected data with correction values, 0 otherwise.

**The following output files are generated:**

- ``fend_object.hdf5``
- ``HiC_data_object.hdf5``
- ``HiC_project_object.hdf5``
- ``HiC_project_object_with distance_parameters.hdf5``
- ``HiC_norm_binning.hdf5`` to be used in the following section.


## 2.

python2.7 /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/HiCtool_compartment_analysis.py \
--action plot_pc \
-i HiCtool_chr6_1mb_PC1.txt \
-c /mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/HiCtool/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1