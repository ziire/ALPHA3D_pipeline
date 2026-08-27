# Automatizing phenomic research: the ALPHA3D pipeline
A toolset for extracting in a fast and efficient way phenotypic traits. Processing automated and manual landmarks. Calculate precision, bias and accuracy of the pipeline. Apply quantitative genetics method to the dataset.

__Article:__ Unpublished

__Authors:__ Irene Zanandrea<sup>1†</sup> and Kjetil Lysne Voje<sup>1</sup>

__Affiliation:__ <sup>1</sup>Evolution and Paleobiology, Natural History Museum, University of Oslo

__Contact:__ <sup>†</sup>irene.zanandrea@nhm.uio.no

__Journal:__ NA

__Year:__ NA  

__Abstract:__ TBA

__Info:__ This repository contains scripts and data used for analyses in the publication.

__Files__ 

_data –_ in this folder there are the landmarks used in the scripts for the analyses.

_scripts –_ this folder contains the scripts used in the analyses. All scripts are commented to be able to reproduce the results from the article.
<ul>
  <li>1_imprecision_log.R runs imprecision (as variance) analysis for the automated pipeline.</li>
  <li>2_bias_log.R runs bias analyses.</li>
  <li>3_variance_ML_vs_AL_ log.R runs imprecision (as variance) between automated and manual landmarking.</li>
  <li>4_inaccuracy_log.R runs inaccuracy analysis (precision plus bias squared) .</li>
  <li>5_cov_matrix_log.R creates phenotypic variance-covariance matrices (P-matrices).</li>
  <li>6_RS_analysis_log.R runs quantitative genetics analyses (Random Skewers).</li>
  <li>7_CS_ANOVA_log.R calculates centroid sizes (CS) and runs nested ANOVAS and Procrustes ANOVA.</li>
  <li>8_meshdist.R runs the meshdist analysis.</li>
  <li>9_Procrustes_ML.R runs Procrustes ANOVA on manual landmarks only.</li>
  <li>create_lists_and_dataframe.R create lists and dataframe used in the above scripts</li>
  <li>functions_accuracy.R contains R functions loaded in all the scripts.</li>
  <li>imprecision_2nd_part_pipeline.R runs imprecision analysis for the second half of the automated pipeline.</li>
  <li>landmarks_traits.R specifies the name of traits used in the above scripts.</li>
  <li>libraries.R contains the libraries used in all the scripts.</li>
  <li>list_nested_anova.R specifies the name of traits used in the ANOVA analyses.</li>
</ul>

