# A/B compartment analysis

This section allows to calculate the principal components (PC) of the Pearson correlation matrix that can be used to delineate A/B compartments in Hi-C data at low resolution (usually 1 mb or 500 kb). It is possible to calculate both PC1 and PC2. Usually, the sign of the eigenvector (PC1) indicates the compartment.

## Table of contents

1. [Calculating the principal component](#1-calculating-the-principal-component)
2. [Plotting the principal component](#2-plotting-the-principal-component)

## 1. Calculating the principal component

HiCtool allows to calculate either the first (typically used) or the second principal component of the Person correlation matrix. In order to do so, the Person correlation matrix has to be calculated first as presented [here](/tutorial/normalization-yaffe-tanay.md#22-normalizing-enrichment-oe-data-and-calculating-the-pearson-correlation-matrix).

This is the code to calculate the first principal component:
```unix
python2.7 /HiCtool-master/scripts/HiCtool_compartment_analysis.py \
--action calculate_pc \
-c /HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1
```

where:

- ``--action``: Action to perform (here ``calculate_pc``).
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--chr``: The chromosome to be used.
- ``--pc``: Which principal component to be returned (either ``PC1`` or ``PC2``).

The output is a txt file with the values of the principal component selected, in this case ``chr6_1000000_PC1.txt`` inside the folder ``./yaffe_tanay_1000000/``.

The parameter ``--chr`` can be used also to pass multiple chromosomes (as a list between square brackets) at once and also multi-processing computation is provided if your machine supports it using the parameters ``-p``:
```unix
chromosomes=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y]

python2.7 /HiCtool-master/scripts/HiCtool_compartment_analysis.py \
--action calculate_pc \
-c /HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr $chromosomes \
--pc PC1
```

## 2. Plotting the principal component

The following code allows to plot the principal component values in a barplot, which allows to visually delineate compartments. The following plot usually is presented together with the Pearson correlation matrix.

```unix
python2.7 /HiCtool-master/scripts/HiCtool_compartment_analysis.py \
--action plot_pc \
-i ./yaffe_tanay_1000000/chr6_1000000_PC1.txt \
-c /HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--pc PC1
```

where:

- ``--action``: Action to perform (here ``plot_pc``).
- ``-i``: Input file with the principal component values calculated above.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--chr``: The chromosome to be used.
- ``--pc``: Which principal component to be returned (the same you have calculated above).

![](/figures/HiCtool_chr6_1mb_PC1.png) 