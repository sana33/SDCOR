# LOF & LoOP

## Two state-of-the-art density-based techniques for outlier detection

In density-based methods, the local density of each object is calculated in a specific way and then is utilized to define the outlier scores. Given an object, the lower its local density compared to that of its neighbors, the more likely it is that the object is an outlier. Density around the points could be calculated by using many techniques, which most of them are distance-based. For example, [Breunig][1] propose a Local Outlier Factor (LOF) that uses the distance values of each object to its nearest neighbors to compute local densities. However, LOF has a drawback which is that the scores obtained through this approach are not globally comparable between all objects in the same dataset or even in different datasets. [Kriegel et al.][2] introduce the Local Outlier Probability (LoOP), which is an enhanced version of LOF. LoOP gives each object a score in the interval [0,1], which is the probability of the object being an outlier and is widely interpretable among various situations.

[1]: Breunig, Markus M., et al. "LOF: identifying density-based local outliers." Proceedings of the 2000 ACM SIGMOD international conference on Management of data. 2000.

[2]: Kriegel, Hans-Peter, et al. "LoOP: local outlier probabilities." Proceedings of the 18th ACM conference on Information and knowledge management. 2009.

