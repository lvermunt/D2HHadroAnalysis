#!/bin/bash

python3 do_training.py -r MLTraining/default_complete.yml -d database/PaperProposal_LcPbPb/database_ml_parameters_LcpK0s3050_2018_binary_fine_final.yml
nice python3 do_application.py -r CreateProcessTTree/default_complete.yml -d database/PaperProposal_LcPbPb/database_ml_parameters_LcpK0s3050_2018_binary_fine_final.yml

