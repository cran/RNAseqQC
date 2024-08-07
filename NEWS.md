# RNAseqQC 0.2.1
* bug fix to not convert factors to numerical before PCA plotting

# RNAseqQC 0.2
* Add the possibility to plot gene expression in plot_pca()
* plot_pca() and plot_pca_scatters() now detect if the variable to color by is
continuous or discrete and chooses an appropriate, colorblind friendly color scale
* Add option to plot rasterised points in plot_pca_scatters()
* It is now possible to plot squared loadings in the plot_loadings() function
* Fix a bug where the distance parameter was ignored in plot_sample_clustering()
* Fix a bug where the annotate_top_n parameter didn't work in plot_loadings() when
the input was a matrix

# RNAseqQC 0.1.4
* Fix of bug in plot_pca() that was introduced in the previous version, causing
the the feature selection to be circumvented.
* Add three dots argument to plot_sample_clustering() to modify the heatmap
* Add documentation to plot_sample_clustering() indicating that it can be used with
an arbitrary SummarizedExperiment object.
* Improved truncation of too large values in plot_chromosome()

# RNAseqQC 0.1.3
* MA plots between replicates better handle missing values
* chromosome heatmaps are not scaled by default anymore
* plot_pca better handles missing values by first mean-imputing NAs and then selecting top-variable features
* added a function to make a matrix of PCA scatter plots to plot each PC against each other
* allow to specify a design during make_dds

# RNAseqQC 0.1.2
* Fixed issue when the results object fed into plot_ma() contains s-values.
* Add GC content gene annotation. 

# RNAseqQC 0.1.1
* Fixed some package dependencies by moving ggsci to Imports and updating the data vignette.

# RNAseqQC 0.1.0
* Initial version
