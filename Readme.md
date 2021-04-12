# SDCOR

## Scalable Density-based Clustering for Local Outlier Detection in Massive-Scale Datasets

#### *Authors:* Sayyed-Ahmad Naghavi-Nozad, Maryam Amir Haeri and Gianluigi Folino

[![arXiv e-print](https://img.shields.io/badge/arXiv-e--print-blue?style=for-the-badge&logo=arXiv&logoColor=white)](https://arxiv.org/pdf/2006.07616.pdf) &emsp; [![Elsevier e-print](https://img.shields.io/badge/Elsevier-e--print-red?style=for-the-badge&logo=Elsevier&logoColor=white)]()

[![license](https://img.shields.io/github/license/DAVFoundation/captain-n3m0.svg?style=flat-square)](https://github.com/sana33/SDCOR/blob/master/LICENSE)

#### Table of Contents

* [Abstract](#abstract)
* [Framework](#framework)
* [Implementations Description](#-implementations-description)
  * [SDCOR _ with visualization - read data from RAM](#-sdcor-_-with-visualization---read-data-from-ram)
    * ["SDCOR Params" panel](#sdcor-params-panel)
    * ["DBSCAN Param Choosing" panel](#dbscan-param-choosing-panel)
    * ["Final Results" panel](#final-results-panel)
    * ["SDCOR various states" panel](#sdcor-various-states-panel)
    * [Main buttons of the GUI](#main-buttons-of-the-gui)
  * [SDCOR _ without visualization - read data from Disk](#-sdcor-_-without-visualization---read-data-from-disk)
    * ["MaxRun" parameter](#maxrun-parameter)
    * ["Progression Status" panel](#progression-status-panel)
    * ["Final Results" panel](#final-results-panel-1)
* [Some Notes](#-some-notes)
* [Citation](#-citation)

#### Abstract:

This paper presents a batch-wise density-based clustering approach for local outlier detection in massive-scale datasets. Unlike the well-known traditional algorithms, which assume that all the data is memory-resident, our proposed method is scalable and processes the input data chunk-by-chunk within the confines of a limited memory buffer. A temporary clustering model is built at the first phase; then, it is gradually updated by analyzing consecutive memory-loads of points. Subsequently, at the end of scalable clustering, the approximate structure of the original clusters is obtained. Finally, by another scan of the entire dataset and using a suitable criterion, an outlying score is assigned to each object called SDCOR (Scalable Density-based Clustering Outlierness Ratio). Evaluations on real-life and synthetic datasets demonstrate that the proposed method has a low linear time complexity and is more effective and efficient compared to best-known conventional density-based methods, which need to load all data into the memory; and also, to some fast distance-based methods, which can perform on data resident in the disk.

#### Framework:

In more detail, the proposed method consists of three major phases. In the first phase, a preliminary random sampling is conducted in order to obtain the main premises on which the algorithm works, i.e., some information on the original clusters and some parameters useful for the progressive clustering. In the second phase, a scalable density-based clustering approach is carried out in order to recognize dense areas on the basis of the currently loaded chunk of data points in the memory. Clusters built incrementally in this way are called miniclusters or subclusters, and they form the temporary clustering model. In more detail, after loading each chunk of data, according to the points already loaded in the memory and those undecided from the previous chunks, and by employing the Mahalanobis distance measure and in respect to the density-based clustering criteria, we update the temporary clustering model, which consists of making some changes to the existing miniclusters or adding new subclusters.

Note that, in the whole scalable clustering procedure, our endeavor aims not to let outliers participate actively in forming and updating any minicluster, and thus, after processing the entire chunks, there will be some objects in the buffer remained undecided. Some of these data are true outliers, while others are inliers, which, due to constraints, have failed to play an influential role in forming a subcluster. Finally, all these undecided points are cleared from the buffer, while only the structural information of the temporary clusters is maintained in the memory. Then, at the last part of scalable clustering, depending on the miniclusters associated with each initial minicluster out of the ``Sampling'' stage, we combine them to obtain the final clusters, which their structure will be approximately the same as of the original clusters.

At last, in the third phase of the proposed approach, w.r.t. the final clustering model gained out of the second phase, once again, we process the entire dataset in chunks to give each object an outlying score, according to the same Mahalanobis distance criterion. The following figure, illustrates the software architecture of the approach.

![Software Architecture](/images/SoftArch.jpg)

## &#10002; Implementations Description

The code is implemented using MATLAB R2016b (version 9.1), and all the experiments are executed on a laptop having a 2.5 GHz Intel Core i5 processor, 6 GB of memory, and a Windows 7 Ultimate (Service Pack 1) operating system. There are two independent versions of implementations, each equipped with a sophisticated GUI:

### &#x1F537; __*SDCOR _ with visualization - read data from RAM*__

![SDCOR with visualizations](/images/SDCOR_RAMversion.png)

In this version, for the ease of visualizations, the input data along with the anomaly labels are totally loaded into memory. Therefore, various kinds of plots could be provided, and moreover, for different steps of the algorithm, there are facilities to plot the data with specific parameters. The details of the GUI are presented as follows:

* ### "SDCOR Params" panel

  * **ChunkSize:** Maximum number of objects in each chunk.
  * **PCvarRatio(%):** PC total variance ratio (in percentage terms).
  * **Alpha:** Membership threshold.
  * **Beta:** Pruning threshold.
  * **SampRate(%):** Random sampling rate (in percentage terms).

  * **Top-n OLs:** Number of top-n outliers for being depicted in the main axes plot.
  * **ScorDSszCoef:** The coefficient value for the obtained outlierness scores to be represented more viewable in the main axes plot.

  * **BlckSzLim:** This parameter is just for expediting the process of the DBSCAN algorithm which is employed in SDCOR. As MATLAB R2016b does not support DBSCAN with a fast built-in C++ function, like for K-means, and more importantly, because MATLAB performs considerably slower than the other programming frameworks in the case of the FOR/WHILE loops, thus we decided to implement our own version of the DBSCAN algorithm.
  
    This version of the DBSCAN implementation might not be so efficient, though it works pretty well on large datasets. DBSCAN has [two ways](https://en.wikipedia.org/wiki/DBSCAN) to be implemented. One is the query-based version which needs to be done in multiple iterations, and hence, is like a poison to MATLAB!; the other one is based on the *Neighbor Graph* gained out of the *n-by-n* distance matrix of the entire data (*n* stands for the cardinality of the input data), which shall be calculated at the first place. We choose the second way, as there is a fast C++-based built-in function in MATLAB, named *pdist2()*, for computing the pairwise distances of the whole objects in data.
  
    Although when the size of the input data goes so high, then the output distance matrix will become too large which even sometimes can not be fit into memory. Moreover, we need not the entire distance matrix to be created at first, and then go for obtaining the *Neighbor Graph*, but we can acquire the distance matrix in small blocks, and then convert each distance block to the corresponding graph block. This could be done by turning each element of the distance block which has a value less than or equal to the *Eps* parameter of DBSCAN to 1, and 0 otherwise.
  
    The blocks are in square shape, and the *BlckSzLim* is the length of the square side. Besides, there is no need for *n* to be divisible by *BlckSzLim*, as our devised algorithm can handle it. Finally, as every element of the distance block is of the "double" type equal to 8 bytes in MATLAB, hence you should consider the *usual* free space of your RAM buffer and then set a reasonable value for this parameter. For example, if the memory free space is equal to 1 GB, then it would be better to consider, e.g., 0.7 GB, for the distance block, which leads to _BlckSzLim = sqrt((0.7×2^30)/8) ≈ 9692_, and leave some space for other operations. The bigger size for _BlckSzLim_, the faster the density-based clustering process will be carried out.
	
* ### "DBSCAN Param Choosing" panel
	
  * #### "Mode" sub-panel
  
	* **Kgraph** Use the _k_-distance graph introduced in the original paper of DBSCAN to locate the optimal parameters.

    * **PSO** Set the PSO evolutionary algorithm for finding the optimal parameters of DBSCAN to operate on the sampled data.
	
    * **Manual** Set the DBSCAN parameters manually to operate on the sampled data.
	
	  _**Note:**_ Please beware that if the structural characteristics of the input data is available, then it would be possible to acquire the optimal parameters of DBSCAN manually. However, such a state is out of our concern here.

  * #### "Initial Params" sub-panel

    * **particleNo:, maxIter:, W:, C1:, C2:, Alpha:** Parameters of the PSO algorithm that you can leave them as default.
	  
	  _**Note:**_  You can even modify the source code and change the current cost function of the PSO algorithm to an arbitrary one.

    * **manuEps:** The manual value for the *Eps* parameter of DBSCAN, set by the user.
    * **manuMnPt:** The manual value for the *MinPts* parameter of DBSCAN, set by the user.

    * **epsCoef:** The coefficient value for the *Eps* parameter to be used while clustering the original distribution; you can leave it as suggested.

  * #### Axes plot

    This plot is for showing the variations of the cost function employed by the PSO algorithm.

  * #### "Make Manu" button

    This button is active when the **Mode** is set to *PSO*. By pressing this button, the optimal parameter values obtained out of PSO algorithm will be set as manual; and thus, in the next run of the proposed method, there will no time spent on finding the optimal values for DBSCAN parameters, to be used for the sampled data.

  * #### "origK:" static text box

    After DBSCAN is applied to the sampled data, the distinct value for the number of original clusters in the input data is attained, which will be displayed in this text box; and will be utilized in the upcoming steps of the proposed method.

* ### "Final Results" panel

  * #### "Accuracy per Chunk" plot
  
    This plot shows the gradual progress of the scalable clustering algorithm in terms of distinguishing outliers in each memory process. The progression is represented w.r.t. various accuracy measures, namely _AUC_, _Precision_, _Recall_ and _F1-Measure_. After processing each chunk of data, we expect the number of inliers in the retained set to become decreased, and on the other side, the number of true outliers become increased. Thus, the mentioned accuracy measures should rise gradually after processing each memory load of points and reach to the perfect condition, iff the parameters are correctly set and also, the input data follows the predefined strong assumptions of the proposed method.
	
	However, after processing the whole chunks, there could be some inliers still sustained in buffer, because of the established membership and density-based clustering restrictions; and moreover, some of the outliers might be absorbed to some created mini-clusters during the scalable clustering. Finally, all the undecided points in memory are cleared from the buffer, while only the structural information of the temporary clusters is maintained in it.
	
	Furthermore, for the density-based anomaly detection algorithms, this axes plot will show the progress as a bar chart which is updated after processing each block of pairwise distances.
	
  * #### "Final AUC:" static text box
  
    This field shows the final AUC obtained out of the "Scoring" phase of the proposed method.
  
  * #### "Time(sec):" static text box
  
    The execution time of SDCOR, regardless of the time spent on visualization matters, is displayed in this field.
  
* ### "SDCOR various states" panel

  * #### Main plot
    
	In this axes plot, different states of the proposed method are illustrated.

  * #### "DispPlot" checkbox
    
	If this checkbox is on, then various states of SDCOR are displayed in the main plot.
	
  * #### "AuxiFig" checkbox

    If this checkbox is on, then each state will be plotted in a specific auxiliary Figure window, instead of the main plot.
	
  * #### "LabelDS" button

    This button colorfully displays the input data including both inliers and outliers.
	
  * #### "SampleDS" button
  
    This button depicts the result of DBSCAN on the sampled data (including noise points) in colors.
	
  * #### "RetSetRed" button
  
    This button illustrates inliers and outliers in different colors, while those points sustained in the retained set after processing the entire chunks and before building the final clustering model, are represented in red color. Besides, the temporary centroids of the temporary clustering model are denoted as magenta square points.
  
  * #### "FinalMeans" button
  
    This button demonstrates both temporary means and final means, represented by solid circles and triangles respectively, with a different color for each final cluster.
  
  * #### "RegenDS" button
  
    This button colorfully demonstrates pruned regenerated data points for every final cluster besides the final means, denoted as dots and triangles respectively.
  
  * #### "ScoredDS" button
  
    This button colorfully illustrates the outcome of "Scoring" phase of the proposed method, based on the local Mahalanobis distance. Considering the final clustering model attained out of the second phase of SDCOR, each point regardless of being an inlier or outlier is denoted as a dot with the color of the closest final cluster, and with a size equal to the corresponding local Mahalanobis distance.
	
	_**Note:**_  The *ScorDSszCoef* input parameter in the *SDCOR Params* panel of the GUI can be employed here to better visualize the scored data points. The bigger this coefficient, the larger the scored points are represented.
  
  * #### "Top-n OLs" button
  
    This button displays the Top-n outliers regarding the acquired outlierness scores out of SDCOR, as red dots. Any other point is represented with a blue dot. The value for _n_ could be modified through the _Top-n OLs_ input parameter in the *SDCOR Params* panel.

* ### Main buttons of the GUI

  * #### "LOAD" button
  
    This button loads the input _n-by-p_ matrix of dataset _X_, along with the _n-by-1_ vector _y_ of outlier labels, all together as a single binary MAT-file. The name of the chosen dataset is depicted in the bold red box, positioned right below the main buttons.
	
	_**Note:**_ In vector _y_ of outlier labels, 0 and 1 stand as the label of inliers and outliers respectively.
	
  * #### "START" button
  
    This button starts the SDCOR algorithm. After pressing this button, any active element of the GUI is disabled, and when the process finishes, all of them will be re-enabled.
  
  * #### "CLEAR AXES" button
  
    This button clears all the axes plots along with the static text fields dedicated to depict a property or the outputs of the selected algorithm.
  
  * #### "RESET BTNS" button
  
    If during the operation, any error happens and the process is halted, or it is stopped intentionally, then all the disabled elements remain inactive. This button activates them again.
  
  * #### "SAVE WORKSPACE" button
  
    This button saves the whole workspace including the output results of the selected algorithm to a single MAT-file, with a suggested file name while saving the workspace.
  
  * #### "LOAD WORKSPACE" button
  
    This button loads a saved workspace from before into the GUI.

### &#x1F536; __*SDCOR _ without visualization - read data from Disk*__

![SDCOR without visualizations](/images/SDCOR_DiskVersion.png)

In this version, the input data is directly read from the disk, and hence, any arbitrary size of dataset could be employed by the proposed method. Therefore, the dedicated facilities for visualizing different steps of the scalable clustering algorithm are removed from the GUI; and instead, some text fields are added to simply demonstrate the algorithm progression and detection accuracy outcomes. The new parts added to the GUI are characterized as follows:

* ### "MaxRun" parameter

   This parameter is added to the *SDCOR Params* panel, and stands for the maximum number of times that SDCOR will be executed on the input data. The final accuracy outcome will be an aggregate value of the total independent runs.
  
* ### "Progression Status" panel

  * #### "Progress by Chunk/Block:" static text box
  
    This field shows the gradual progress of the scalable clustering algorithm or the chosen density-based anomaly algorithm, in terms of successive chunks or blocks of data, respectively.
	
  * #### "Temp AUC:" static text box
  
    This field displays the AUC outcome obtained through the last run of the proposed method, out of the _MaxRun_ number of execution times.
	
  * #### "Temp Time(sec):" static text box

    This field displays the execution time (in seconds) of the last run of the proposed method, out of the _MaxRun_ number of execution times.
	
  * #### "Total Runs:" static text box
  
    The total number of runs of SDCOR per the allowed maximum number of execution times (_MaxRun_) is depicted in this field.
	
* ### "Final Results" panel

  * #### "Final AUC:" static text box
  
    This field shows the average AUC outcome attained through the _MaxRun_ number of execution times of SDCOR.
	
  * #### "AUC std:" static text box

    This field shows the standard deviation of AUC outcomes gained out of the _MaxRun_ number of SDCOR runs.
	
  * #### "Time(sec):" static text box
  
    The average execution time (in seconds) acquired through the _MaxRun_ number of SDCOR runs.
	
## &#x1F34F; Some Notes

* #### For the better experience of working with MATLAB GUI, right-click on the FIG-file and select "Open in GUIDE", instead of double-clicking on it.

* #### There are some commented scripts (each followed by _"@-- debugging script --@"_) in the code that could be uncommented if necessary, in the case of tracking the functionality of the proposed method.

* #### There are two training videos for both [RAM version](/videos/SDCOR_RAMversion.wmv) and [Disk version](/videos/SDCOR_DiskVersion.wmv), which each of them demonstrates a test run of the corresponding implementation.

* #### If the process of a GUI gets so long and intolerable, or if for any kind of reason you want to cease the operation, you can click in the _Command Window_ of MATLAB and press _Ctrl+C_ or _Ctrl+Break_.

## &#x1F4D8; Citation

If you are interested in the idea or you are using this code for your research, please cite our paper as:

```
@article{naghavi2020sdcor,
  title={SDCOR: Scalable Density-based Clustering for Local Outlier Detection in Massive-Scale Datasets},
  author={Naghavi-Nozad, Sayyed-Ahmad and Haeri, Maryam Amir and Folino, Gianluigi},
  journal={arXiv preprint arXiv:2006.07616},
  year={2020}
}
```

Thanks a lot ... :pray:&#x1F49C;