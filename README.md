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
  <li>single_mode.R runs the analyses for the single-mode models section of the article.</li>
  <li>stats_single_mode.R runs regression analyses and plotting on the single-mode data generated in single_mode.R.</li>
  <li>mode_shift.R runs the analyses pluss summary statistics for the mode-shift models section of the article.</li>
  <li>stats_mode_shift.R runs regression analyses and plotting on the mode-shift data generated in mode_shift.R.</li>
  <li>shift_OU_adequacy.R is a script run on a HPC cluster for the OU part of adequacy testing in the mode-shift part of the article.</li>
  <li>delta_aicc_gap.R runs analyses and plotting for the delta AICc gap part of the article.</li>
  <li>functions.R contains R functions loaded in all the above scripts.</li>
</ul>

