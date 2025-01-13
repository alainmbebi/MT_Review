# MT_Review

R functions for the comparative analysis of: Advances in multi-trait genomic prediction approaches: Classification, comparative analysis, and perspectives, as described in [( Mbebi A. et al. 2025)]( https://doi.org/10.1).

1.    The folder code contains the following R scripts:
  
*    Genotype_BLUP_Metabolites.R which computes the BLUP for each metabolite trait
*    Focal_Trait_Heritability.R which computes the heritability for focal traits
*    Ridge_estim.R compute the Ridge estimate
*    Qr_Reorthogonalization.R for re-othogonalization of QR decomposition
*    ST_GBLUP_CVA.R which implements the single-trait GBLUP for the scenario where the model is trained on Indica to predict Indica
*    MT_SVD_CVA.R which implements the MT SVD for the scenario where the model is trained on Indica to predict Indica
*    MT_PLS_CVA.R which implements the MT PLS for the scenario where the model is trained on Indica to predict Indica
*    MT_MOR_CVA_Indica_Indica.R which implements the MT MOR for the scenario where the model is trained on Indica to predict Indica
*    MT_DL_CVA.R which implements the MT deep learning (MT-DL) for the scenario where the model is trained on Indica to predict Indica
*    MT_BMORS_CVA.R which implements the MT BMORS for the scenario where the model is trained on Indica to predict Indica

  
3.    The folder data contains all the data used in the article
4.    Notes
*  Although the codes here were tested on Ubuntu 22.04.5 LTS using R (version 4.2.2), they can run under any Linux or Windows OS distributions, as long as all the required packages are compatible with the desired R version.
*  The following abbreviations are used, ST: Single-trait MT: Multi-trait, DL: Deep learning, SVD: Singular value decomposition, GBLUP: genomic best linear unbiased prediction, PLS: partial least square regression, MOR: Multi output regression, BMORS: Bayesian multi-output regressor stacking, CVA: Cross validation A (see manuscript for details).

5.  Licence: GPL-3
