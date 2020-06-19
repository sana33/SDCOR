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

### &#x1F537; __*SDCOR _ with visualization - read data from RAM*__:

![SDCOR with visualizations](/images/SDCOR_RAMversion.png)

In this version, for the ease of visualizations, the input data along with the anomaly labels are totally loaded into memory. Therefore, various kinds of plots could be provided and moreover, for different steps of the algorithm, there are facilities to plot the data with specific parameters. The details of the GUI are represented as follows:

* ### "SDCOR Params" panel

  * **ChunkSize:** Number of objects in each chunk.
  * **PCvarRatio(%):** PC total variance ratio (in percentage terms).
  * **Alpha:** Membership threshold.
  * **Beta:** Pruning threshold.
  * **SampRate(%):** Random sampling rate (in percentage terms).

  * **Top-n OLs:** Number of top-n outliers for being depicted in the main axes plot.
  * **ScorDSszCoef:** The coefficient value for the obtained outlierness scores to be represented more viewable in the main axes plot.

  * **BlckSzLim:** This parameter is just for expediting the process of the DBSCAN algorithm which is employed in SDCOR. As MATLAB 9 does not support the DBSCAN algorithm with a fast built-in C++ function, like K-means; and more importantly, because MATLAB is seriously slow in loops (like 'for' and 'while' loops), thus we decided to implement DBSCAN with a code of our own.
  
    This version of DBSCAN implementation might not be so efficient, though it works pretty well on large datasets. DBSCAN has [two ways](https://en.wikipedia.org/wiki/DBSCAN) to be implemented. One is the query-based version which needs to be done in multiple iterations, and hence, is like a poison to MATLAB!; and the other one is based on the *Neighbor Graph* which is gained out of the *n-by-n* distance matrix of the entire data (*n* stands for the cardinality of the input data), which shall be calculated at the first place. We choose the second way, as there is a fast C++-based built-in function in MATLAB, named *pdist2()*, for computing the pairwise distances of the whole objects in data.
  
    Although when the size of the input data goes so high, then the output distance matrix will become too large which even sometimes can not be fit into memory. Moreover, we need not the entire distance matrix to be created at first, and then go for obtaining the *Neighbor Graph*; but we can acquire the distance matrix in small blocks, and then convert each distance block to the corresponding graph block. This could be done by changing each element of the distance block, which has a distance value less than or equal to the *Eps* parameter of DBSCAN, to 1, and otherwise to 0.
  
    The blocks are in square shape, and the *BlckSzLim* is the length of the square side. Besides, there is no need for *n* to be divisible by *BlckSzLim*, as our devised algorithm can handle it. Finally, as each element of the distance block is of double type, which is equal to 8 bytes in MATLAB; hence, you should consider the *usual* free space of your RAM buffer and then set a reasonable value for this parameter. For example, if the free space in memory is equal to 1 GB, then it would be better to consider e.g. 0.7 GB for the distance block, which leads to _BlckSzLim = sqrt((0.7×2^30)/8) ≈ 9692_, and leave some space for other operations. The bigger size for _BlckSzLim_, the faster the density-based clustering process will be carried out.
	
	Furthermore, for the materialization matrix required by LOF and LoOP, this parameter is used in a slightly different manner.
	
* ### "DBSCAN Param Choosing" panel
	
  * #### "Mode" sub-panel

    * **PSO** Set PSO evolutionary algorithm for finding the optimal parameters of DBSCAN algorithm to operate on the sampled data. Note that in this mode, the spent time on locating the optimal parameters is taken into account for the total runtime of SDCOR.
    * **Manual** Set the DBSCAN parameters manually to operate on the sampled data.
	
	  _**Note:**_ Beware that if the structural characteristics of the utilized input data is available, then it would be possible to acquire the optimal parameters of DBSCAN manually.

  * #### "Initial Params" sub-panel

    * **dimCoef:, particleNo:, maxIter:, W:, C1:, C2:, Alpha:** Parameters of PSO algorithm, which you can leave them as default.
	  
	  _**Note:**_  You can even modify the source code and change the current cost function of PSO algorithm to a better one.

    * **manuEps:** The manual value for the *Eps* parameter of DBSCAN, set by the user.
    * **manuMnPt:** The manual value for the *MinPts* parameter of DBSCAN, set by the user.

    * **epsCoef:** The coefficient value for the *Eps* parameter to be used while clustering the original distribution. You can leave it as suggested.
    * **MinPtsCoef:** The coefficient value for the *MinPts* parameter to be used while clustering the original distribution. You can leave it as suggested.

  * #### Axes plot

    This plot is for showing the variations of the cost function employed by the PSO algorithm.

  * #### "Make Manu" button

    This button is active when the **Mode** is set to *PSO*. By pressing this button, the optimal parameter values obtained out of PSO algorithm will be set as manual; and thus, in the next run of the proposed method, there will no time spent on finding the optimal values for DBSCAN parameters, to be used for the sampled data.

  * #### "origK:" static text box

    After DBSCAN is applied to the sampled data, the distinct value for the number of original clusters in the input data is attained, which will be displayed in this text box; and will be utilized in the upcoming steps of the proposed method.

* ### "Final Results" panel

  * #### "Accuracy per Chunk" plot
  
    This plot shows the gradual progress of the scalable clustering algorithm in terms of distinguishing outliers in each memory process. The progression is represented w.r.t. various accuracy measures, namely _AUC_, _Precision_, _Recall_ and _F1-Measure_. After processing each chunk of data, we expect the number of inliers in the retained set to become decreased, and on the other side, the number of true outliers become increased. Thus, the mentioned accuracy measures should rise gradually after processing each memory load of points and reach to the perfect condition, iff the parameters are correctly set and also, the input data follows the predefined strong assumptions of the proposed method. However, after processing the whole chunks, there could be some inliers still sustained in buffer, because of the established membership and density-based clustering restrictions; and moreover, some of the outliers might be absorbed to some created mini-clusters during the scalable clustering. Finally, all the undecided points in memory are cleared from the buffer, while only the structural information of the temporary clusters is maintained in it.
	
	For the density-based anomaly detection algorithms, this axes plot will show the progress as a bar chart which is updated after processing each block of pairwise distances.
	
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
	
* ### "LOF & LoOP Params" panel

  * #### "LOF" checkbox
    
    If this checkbox is on, then the LOF algorithm will be applied to the input data, considering the corresponding input parameters. By clicking this checkbox, all of the other elements of the GUI associated with other competing methods will be disabled.
    
  * #### "LoOP" checkbox

    If this checkbox is on, then the LoOP algorithm will be applied to the input data, considering the corresponding input parameters. By clicking this checkbox, all of the other elements of the GUI associated with other competing methods will be disabled.
	
  * **MinPts Intv.:** This parameter is modified w.r.t. the selected density-based outlier algorithm, namely LOF or LoOP. The default established values are set as suggested. For LOF, it is an interval, and for LoOP, it becomes a scalar value.

  * **k StepLength:** Step-length for the MinPts interval utilized by LOF.

  * **Lambda:** Lambda parameter of the LoOP algorithm.

* ### Main buttons of the GUI

  * #### "LOAD" button
  
    This button loads the input _n-by-p_ matrix of dataset _X_, along with the _n-by-1_ vector _y_ of outlier labels, all together as a single binary MAT-file. The name of the chosen dataset is depicted in the bold red box, positioned right below the main buttons.
	
  * #### "START" button
  
    This button starts the operations of the selected algorithm. If _LOF_ or _LoOP_ checkboxes are not checked, then the SDCOR algorithm will be commenced. After pressing this button, any active element of the GUI is disabled, and when the process finishes, all of them will be re-enabled.
  
  * #### "CLEAR AXES" button
  
    This button clears all the axes plots along with the static text fields dedicated to depict a property or the outputs of the selected algorithm.
  
  * #### "RESET BTNS" button
  
    If during the operation, any error happens and the process is halted, or it is stopped intentionally, then all the disabled elements remain inactive. This button activates them again.
  
  * #### "SAVE WORKSPACE" button
  
    This button saves the whole workspace including the output results of the selected algorithm to a single MAT-file, with a suggested file name while saving the workspace.
  
  * #### "LOAD WORKSPACE" button
  
    This button loads a saved workspace from before into the GUI.

### &#x1F537; __*SDCOR _ without visualization - read data from Disk*__: 

![SDCOR without visualizations](/images/SDCOR_DiskVersion.png)

In this version, the input data is directly read from the disk, and hence, any arbitrary size of dataset could be employed by the proposed method. Therefore, the dedicated facilities for visualizing different steps of the scalable clustering algorithm are removed from the GUI; and instead, some text fields are added to simply demonstrate the algorithm progression and detection accuracy outcomes. The new parts added to the GUI are characterized as follows:

* **MaxRun:** This field is added to the *SDCOR Params* panel, and stands for the maximum number of times that SDCOR will be executed on the input data. The final accuracy outcome will be an aggregate value of the total independent runs.
  
* ### "Progression Status" panel

  * #### "Progress by Chunk/Block:" static text box
  
    This field shows the gradual progress of the scalable clustering algorithm or the density-based anomaly algorithms, in terms of successive chunks or blocks of data, respectively.
	
  * #### "Temp AUC:" static text box
  
    This field displays the AUC outcome obtained through the last run of the proposed method, out of the _MaxRun_ number of execution times.
	
  * #### "Temp Time(sec):" static text box

    This field displays the execution time (in seconds) of the last run of the proposed method, out of the _MaxRun_ number of execution times.
	
  * #### "Total Runs:" static text box
  
    The total number of runs of SDCOR per the maximum number of execution times (_MaxRun_) is depicted in this field.
	
* ### "Final Results" panel

  * #### "Final AUC:" static text box
  
    This field shows the average AUC outcome attained through the _MaxRun_ number of execution times of SDCOR.
	
  * #### "AUC std:" static text box

    This field shows the standard deviation of AUC outcomes gained out of the _MaxRun_ number of SDCOR runs.
	
  * #### "Time(sec):" static text box
  
    The average execution time (in seconds) acquired through the _MaxRun_ number of SDCOR runs.
	
### &#x1F34F; **_Some Notes:_**

* #### For the better experience of working with the GUI files, right-click the FIG-file in MATLAB and select "Open in GUIDE".

* #### There are some commented scripts (each followed by _"@-- debugging script --@"_) in the code that could be uncommented if necessary, in the case of tracking the performance of the proposed method.

* #### There are two videos for both [RAM version](/videos/SDCOR_RAMversion.wmv) and [Disk version](/videos/SDCOR_DiskVersion.wmv), which each of them demonstrates a test run of the corresponding implementation.

* #### If the process of a GUI gets so long and intolerable, or if for any kind of reason you want to cease the operation, you can click in the _Command Window_ of MATLAB and press _Ctrl+C_ or _Ctrl+Break_.

