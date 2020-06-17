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

The code is implemented using MATLAB 9, and all the experiments are executed on a laptop having a 2.5 GHz Intel Core i5 processor and 6 GB of memory. The implementations of the two state-of-the-art density-based outlier detection methods, namely LOF and LoOP, are included in this code too.

There are two independent versions of implementations, each equipped with a sophisticated GUI:

* __SDCOR _ with visualization - read data from RAM__:

![SDCOR with visualizations](/images/SDCOR_RAMversion.png)
	
This one is provided with various kinds of plots for the ease of visualizations, as for different steps of the algorithm, there are facilities to plot the data with specific parameters. The details of the GUI are represented as follows:

* ### "SDCOR Params" panel:

 * **ChunkSize:** Number of objects in each chunk.
**PCvarRatio(%):** PC total variance ratio (in percentage terms).
**Alpha:** Membership threshold.
**Beta:** Pruning threshold.
**SampRate(%):** Random sampling rate (in percentage terms).

**Top-n OLs:** Number of top-n outliers for being shown in the main axes plot.
**ScorDSszCoef:** The coefficient for the obtained outlierness scores to be represented more viewable in the main axes plot.

**BlckSzLim:** This parameter is just for expediting the process of the DBSCAN algorithm which is employed in SDCOR. As MATLAB 9 does not support the DBSCAN algorithm with a fast built-in C++ function, like K-means; and more importantly, because MATLAB is seriously slow in loops (like 'for' and 'while' loops), thus we decided to implement DBSCAN with a code of our own.
However, this version might not be so efficient, it works pretty well on large datasets. DBSCAN has [two ways](https://en.wikipedia.org/wiki/DBSCAN) to be implemented. One is the query-based version which needs to be done in multiple iterations, and is like a poison to MATLAB!; and the other one is based on the *Neighbor Graph* which is gained out of the *n-by-n* distance matrix of the entire data (*n* stands for the cardinality of the input data), which shall be calculated at the first place. We choose the second way, as there is a fast C++-based built-in function in MATLAB, named pdist2(), for computing the pairwise distances of the whole objects in data.
Although when the size of the input data goes so high, then the output distance matrix will become too large which even sometimes can not be fit into memory. Moreover, we need not the entire distance matrix to be created at first, and then go for obtaining the Neighbor Graph; but we can acquire the distance matrix in small blocks, and then convert each block to the corresponding block of the Neighbor Graph. This could be done by changing each element of the distance block, which has a distance value less than or equal to *Eps* parameter of the DBSCAN, to 1, and to 0 otherwise.
The blocks are in square shape, and the *BlckSzLim* is the length of the square side. Besides, there is no need for *n* to be divisible by *BlckSzLim*, as our devised algorithm can handle it. Finally, as each element of the distance block is of the double type, which is equal to 8 bytes in MATLAB; hence, you should consider the *usual* free space of your RAM buffer and then set a reasonable value for this parameter. For example, if the free space in memory is equal to 1 GB, then it would be better to consider e.g. 0.7 GB for the distance block, which leads to _BlckSzLim = √(0.7×2^30)/8 ≈ 9692_, and leave some space for other operations. The bigger size for the _BlckSzLim_, the faster the density-based clustering process will be carries out.
	
* ### "DBSCAN Param Choosing" panel
	
#### "Mode" sub-panel

**PSO** Set PSO evolutionary algorithm for finding the optimal parameters of DBSCAN algorithm to operate on the sampled data.
**Manual** Set the DBSCAN parameters manually to operate on the sampled data.

#### "Initial Params" sub-panel

**dimCoef:, particleNo:, maxIter:, W:, C1:, C2:, Alpha:** Parameters of PSO algorithm, which you can leave them as default.

**manuEps:** The manual value for the *Eps* parameter of DBSCAN, set by the user.
**manuMnPt:** The manual value for the *MinPts* parameter of DBSCAN, set by the user.

**epsCoef:** The coefficient value for the *Eps* parameter to be used while clustering the original distribution. You can leave it as suggested by the author.
**MinPtsCoef:** The coefficient value for the *MinPts* parameter to be used while clustering the original distribution. You can leave it as suggested by the author.

#### Axes Plot

This plot is for showing the variations of the cost function employed by the PSO algorithm.

#### "Make Manu" button

This button is active when the **Mode** is set to *PSO*. By pressing this button, the optimal parameter values obtained out of PSO algorithm will be set as manual; and thus, in the next run of the proposed method, there will no time spent on finding the optimal values for DBSCAN parameters, to be used for the sampled data.

#### "origK" static text box

After DBSCAN is applied to the sampled data, the distinct value for the number of the original clusters in the input data is attained, which will be displayed in this text box; and will be utilized in the upcoming steps of the proposed method.










* __SDCOR _ without visualization - read data from Disk__: 

![SDCOR without visualizations](/images/SDCOR_DiskVersion.png)

This one is provided with va...







