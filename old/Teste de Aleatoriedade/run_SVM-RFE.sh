#!/bin/bash
# Run Double-Cross-Validation and Rankings-Validation using SVM-RFE and SVM
R -f script_secretoma_SVM_RFE_DoubleCrossValidation.R
R -f script_rank_validation_svmrfe_test_train.R

