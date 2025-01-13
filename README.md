# MT_Review

R functions for the comparative analysis of: Advances in multi-trait genomic prediction approaches: Classification, comparative analysis, and perspectives, as described in ( Mbebi A. et al. 2025).

1.    The folder code contains the following R scripts:

    L21_featselect.R which uses the model proposed in (Nie et al. 2010) to implement GS under the equation-norm regularized multivariate regression

    CV_L21_featselect.R which selects the tuning parameter for L21_featselect.R using K-folds cross-validation

    Ridge_estim.R compute the Ridge estimate

    CV_Ridge.R selects the tuning parameter for Ridge_estim.R using K-folds cross-validation

    L21_joint_estim.R performs GS using the equation-norm regularized multivariate regression that jointly estimates the regression coefficients and precision matrix

    CV_L21_joint_estim.R selects the tuning parameters for L21_joint_estim.R using K-folds cross-validation

    MOR.R runs multiple output regression (He et al. 2016)

    CV_MOR.R selects the tuning parameters for MOR.R using K-folds CV

    2.    The folder data contains all the data used in the article.

    Notes

    Although the codes here were tested on Ubuntu 22.04.5 LTS using R (version 4.2.2), they can run under any Linux or Windows OS distributions, as long as all the required packages are compatible with the desired R version.


    3.    Licence: GPL-3
