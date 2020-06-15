# SDCOR
 SDCOR: Scalable Density-based Clustering for Local Outlier Detection in Massive-Scale Datasets

This paper presents a batch-wise density-based clustering approach for local outlier detection in massive-scale datasets. Differently from well-known traditional algorithms, which assume that all the data is memory-resident, our proposed method is scalable and processes the data chunk-by-chunk within the confines of a limited memory buffer. At first, a temporary clustering model is built, then it is incrementally updated by analyzing consecutive memory loads of points. Ultimately, the proposed algorithm will give an outlying score to each object, which is named SDCOR (Scalable Density-based Clustering Outlierness Ratio). Evaluations on real-life and synthetic datasets demonstrate that the proposed method has a low linear time complexity and is more effective and efficient compared to best-known conventional density-based methods, which need to load all the data into memory; and also some fast distance-based methods which can perform on the data resident in the disk.

## Implementations Description

The code is implemented using MATLAB 9, and all the experiments are executed on a laptop having a 2.5 GHz Intel Core i5 processor and 6 GB of memory.

There are two independent versions of implementations, each equipped with a sophisticated GUI:

* __SDCOR _ with visualization - read data from RAM__:

This one is provided with various axes plots for the ease of visualizations, as for different steps of the algorithm, there are facilities to plot the data with specific parameters.

* __SDCOR _ without visualization - read data from Disk__: 

This one is provided with va...




