#!/bin/bash
# Run Double-Cross-Validation and Rankings-Validation using NSC
R -f script_secretoma_NSC_DoubleCrossValidation.R
R -f script_rank_validation_nsc_test_train.R

