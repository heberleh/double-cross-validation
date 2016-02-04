#!/bin/bash

# ==README==
# These commands run the ranking and validation scripts
# use # to comment and ignore any script
# remove all # from script lines only to run all the scripts.
# Using command line instead of IDE, like RStudio, will increase performance (less time to finish).
# Try it using Ubuntu. Make sure to have installed all the library() necessary.

# Required R libs:
# MASS, pamr, "parallel", ibb, rpart, class, e1071, dismo
# require(foreach), require(doSNOW)

# Instead of changing each input file name inside the scripts
# copy your files to dataset/current and rename them as follow:
# 	> input.txt
# 	> input_beta-binomial.txt

# The results will be stored with generic names as well.
# The will be avaiable at:
#        > ./results/beta-binomial/
#        > ./results/nsc/
#        > ./results/svm-rfe/
# Check the "modified date" of the files to always be sure that
# those results are the ones you want.
# Copy or move the files to a new folder, because each time
# the scripts are executed, the files in /results folder 
# will be replaced by new ones.


# Run Beta-Binomial test
R -f ./src/script_Beta-Binomial.R #>/dev/null

# Run Double-Cross-Validation and Rankings-Validation using NSC
R -f ./src/script_NSC_DoubleCrossValidation.R #>/dev/null
R -f ./src/script_rank_validation_nsc_test_train.R #>/dev/null

# Run Double-Cross-Validation and Rankings-Validation using SVM-RFE and SVM
#R -f ./src/script_SVM_RFE_DoubleCrossValidation.R #>/dev/null
#R -f ./src/script_rank_validation_svmrfe_test_train.R #>/dev/null
