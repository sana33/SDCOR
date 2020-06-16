# SDCOR
## Scalable Density-based Clustering for Local Outlier Detection in Massive-Scale Datasets

#### Abstract:

This paper presents a batch-wise density-based clustering approach for local outlier detection in massive-scale datasets. Differently from well-known traditional algorithms, which assume that all the data is memory-resident, our proposed method is scalable and processes the data chunk-by-chunk within the confines of a limited memory buffer. At first, a temporary clustering model is built, then it is incrementally updated by analyzing consecutive memory loads of points. Ultimately, the proposed algorithm will give an outlying score to each object, which is named SDCOR (Scalable Density-based Clustering Outlierness Ratio). Evaluations on real-life and synthetic datasets demonstrate that the proposed method has a low linear time complexity and is more effective and efficient compared to best-known conventional density-based methods, which need to load all the data into memory; and also some fast distance-based methods which can perform on the data resident in the disk.

#### Framework:

In more detail, the proposed method consists of three major phases. In the first phase, a preliminary random sampling is conducted in order to obtain the main premises on which the algorithm works, i.e. some information on the original clusters and some parameters useful for the incremental clustering. In the second phase, a scalable density-based clustering algorithm is carried out in order to recognize dense areas, on the basis of currently loaded chunk of data points in memory. Clusters built incrementally in this phase, are called mini-clusters or sub-clusters, and they form the temporary clustering model. After loading each chunk of data, according to the points already loaded in memory and those undecided from previous chunks; and by employing the Mahalanobis distance measure and the density-based clustering criteria; we update the temporary clustering model, which consists of making some changes to existing mini-clusters or adding new sub-clusters.

Note that, in the whole scalable clustering procedure, our endeavor is aimed to not let outliers participate actively in forming and updating any mini-cluster, and thus, after processing the entire chunks, there will be some objects in buffer remained undecided. Some of these data are true outliers, while some others are inliers, which, due to constraints, have failed to play any effective role in forming a sub-cluster. Finally, all these undecided points are cleared from the buffer, while only the structural information of the temporary clusters is maintained in memory. Then, at the last part of the scalable clustering algorithm, we utilize another clustering-based approach to combine the mini-clusters and obtain the final clusters, which their structure will be approximately the same as of the original clusters.

At last, in the third phase of the proposed approach, w.r.t. the final clustering model gained out of the second phase, once again, we process the entire dataset in chunks, to give each object an outlying score, according to the same Mahalanobis distance criterion. The following figure, illustrates the software architecture of the approach.

![Software Architecture](/images/SoftArch.jpg)

## Implementations Description

The code is implemented using MATLAB 9, and all the experiments are executed on a laptop having a 2.5 GHz Intel Core i5 processor and 6 GB of memory. The LOF

There are two independent versions of implementations, each equipped with a sophisticated GUI:

* __SDCOR _ with visualization - read data from RAM__:

	![SDCOR with visualizations](/images/SDCOR_RAMversion.png)
	
	This one is provided with various kinds of plots for the ease of visualizations, as for different steps of the algorithm, there are facilities to plot the data with specific parameters. The details of the GUI is represented as follows:
			
	

	

* __SDCOR _ without visualization - read data from Disk__: 

![SDCOR without visualizations](/images/SDCOR_DiskVersion.png)

	This one is provided with va...
	
	
					



