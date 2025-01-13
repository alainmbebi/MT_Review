# MT_Review

R functions for the comparative analysis of: Advances in multi-trait genomic prediction approaches: Classification, comparative analysis, and perspectives, as described in [( Mbebi A. et al. 2025)]( https://doi.org/10.1).

1.    The folder code contains the following R scripts:

    Genotype_BLUP_Metabolites.R which uses the replicates to compute the BLUP for each metabolite trait

    Focal_Trait_Heritability.R which compute the heritability for focal traits

    Ridge_estim.R compute the Ridge estimate

    Qr_Reorthogonalization.R for re-othogonalization of QR decomposition
    
    Plink_steps.PNG contains the preprocessing 

    CV_L21_joint_estim.R selects the tuning parameters for L21_joint_estim.R using K-folds cross-validation

    MOR.R runs multiple output regression (He et al. 2016)

    CV_MOR.R selects the tuning parameters for MOR.R using K-folds CV

    2.    The folder data contains all the data used in the article.

    Notes

    Although the codes here were tested on Ubuntu 22.04.5 LTS using R (version 4.2.2), they can run under any Linux or Windows OS distributions, as long as all the required packages are compatible with the desired R version.


    3.    Licence: GPL-3
