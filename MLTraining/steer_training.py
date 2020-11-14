#############################################################################
##  © Copyright CERN 2018. All rights not expressly granted are reserved.  ##
##                 Author: Gian.Michele.Innocenti@cern.ch                  ##
## This program is free software: you can redistribute it and/or modify it ##
##  under the terms of the GNU General Public License as published by the  ##
## Free Software Foundation, either version 3 of the License, or (at your  ##
## option) any later version. This program is distributed in the hope that ##
##  it will be useful, but WITHOUT ANY WARRANTY; without even the implied  ##
##     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    ##
##           See the GNU General Public License for more details.          ##
##    You should have received a copy of the GNU General Public License    ##
##   along with this program. if not, see <https://www.gnu.org/licenses/>. ##
#############################################################################

"""
main script for doing data processing, machine learning and analysis
"""
import argparse
import sys
from os.path import exists

import warnings
import yaml
from pkg_resources import resource_stream
# To set batch mode immediately
from ROOT import gROOT  # pylint: disable=import-error, no-name-in-module

from machine_learning_hep.logger import configure_logger, get_logger
from machine_learning_hep.utilities import checkdir
from machine_learning_hep.utilities import checkmakedir

#from MLTraining.optimiser import Optimiser
from MLTraining.optimiser_hipe4ml import Optimiserhipe4ml
from MLTraining.optimiser_hipe4mltree import Optimiserhipe4mltree

warnings.simplefilter(action='ignore', category=FutureWarning)

try:
    import logging
    import absl.logging
    logging.root.removeHandler(absl.logging._absl_handler) # pylint: disable=protected-access
    absl.logging._warn_preinit_stderr = False # pylint: disable=protected-access
except Exception as e: # pylint: disable=broad-except
    print("##############################")
    print("Failed to fix absl logging bug", e)
    print("##############################")


def do_training(data_config: dict, data_param: dict, data_model: dict): # pylint: disable=too-many-locals, too-many-statements, too-many-branches

    # Disable any graphical stuff. No TCanvases opened and shown by default
    gROOT.SetBatch(True)

    logger = get_logger()
    logger.info("Do training chain")

    # If we are here we are interested in the very first key in the parameters database
    for k in data_param.keys():
        case = k
        break

    doml = data_config["ml_study"]["activate"]
    domloption = data_config["ml_study"]["dohipe4ml"]
    domlprefilterstep = data_config["ml_study"]["domlprefilterstep"]
    docorrelation = data_config["ml_study"]['docorrelation']
    dotraining = data_config["ml_study"]['dotraining']
    dotesting = data_config["ml_study"]['dotesting']
    doapplytodatamc = data_config["ml_study"]['doapplytodatamc']
    docrossvalidation = data_config["ml_study"]['docrossvalidation']
    dolearningcurve = data_config["ml_study"]['dolearningcurve']
    doroc = data_config["ml_study"]['doroc']
    doroctraintest = data_config["ml_study"]['doroctraintest']
    doboundary = data_config["ml_study"]['doboundary']
    doimportance = data_config["ml_study"]['doimportance']
    dogridsearch = data_config["ml_study"]['dogridsearch']
    doefficiencyml = data_config["ml_study"]['doefficiency']
    dosignifopt = data_config["ml_study"]['dosignifopt']
    doscancuts = data_config["ml_study"]["doscancuts"]
    doplotdistr = data_config["ml_study"]["doplotdistr"]
    domlplots = data_config["ml_study"].get("domlplots", False)

    typean = data_config["analysis"]["type"]

    binminarray = data_param[case]["ml"]["binmin"]
    binmaxarray = data_param[case]["ml"]["binmax"]
    raahp = data_param[case]["ml"]["opt"]["raahp"]
    mltype = data_param[case]["ml"]["mltype"]
    training_vars = data_param[case]["variables"]["var_training"]
    if not isinstance(training_vars[0], list):
        training_vars = [training_vars for _ in range(len(binminarray))]

    domlprefilterstep_indb = data_param[case].get("doml_asprefilter", None)
    if domlprefilterstep is True and domlprefilterstep_indb is not None:
        data_param[case]["doml_asprefilter"] = True
        domlprefilterstep_indb = True
    if domlprefilterstep is False and domlprefilterstep_indb is not None:
        data_param[case]["doml_asprefilter"] = False
        domlprefilterstep_indb = False
    if domlprefilterstep_indb is None:
        data_param[case]["doml_asprefilter"] = None
        domlprefilterstep_indb = None
        if domlprefilterstep is True:
            logger.warning("Mismatch between domlprefilter step in config and DBs")

    mlout = data_param[case]["ml"]["mlout"]
    mlplot = data_param[case]["ml"]["mlplot"]
    if domlprefilterstep_indb is True:
        data_param[case]["ml"]["mlout"] = mlout + "/prefilter"
        data_param[case]["ml"]["mlplot"] = mlplot + "/prefilter"
        mlout = mlout + "/prefilter"
        mlplot = mlplot + "/prefilter"
    if domlprefilterstep_indb is False:
        data_param[case]["ml"]["mlout"] = mlout + "/analysis"
        data_param[case]["ml"]["mlplot"] = mlplot + "/analysis"
        mlout = mlout + "/analysis"
        mlplot = mlplot + "/analysis"

    opti_hyperpar_hipe4ml = data_param[case]["hipe4ml"]["hyper_par_opt"]["do_hyp_opt"]
    hipe4ml_hyper_pars = data_param[case]["hipe4ml"]["hipe4ml_hyper_pars"]

    #creating folder if not present
    counter = 0
    if doml is True and dotraining is True:
        counter = counter + checkdir(mlout)
        counter = counter + checkdir(mlplot)

    if counter < 0:
        sys.exit()

    # check and create directories
    if doml is True and dotraining is True:
        checkmakedir(mlout)
        checkmakedir(mlplot)

#    if doml is True and domloption == 1:
#        index = 0
#        for binmin, binmax in zip(binminarray, binmaxarray):
#            myopt = Optimiser(data_param[case], case, typean,
#                              data_model[mltype], binmin, binmax,
#                              raahp[index], training_vars[index])
#            if docorrelation is True:
#                myopt.do_corr()
#            if dotraining is True:
#                myopt.do_train()
#            if dotesting is True:
#                myopt.do_test()
#            if doapplytodatamc is True:
#                myopt.do_apply()
#            if docrossvalidation is True:
#                myopt.do_crossval()
#            if dolearningcurve is True:
#                myopt.do_learningcurve()
#            if doroc is True:
#                myopt.do_roc()
#            if doroctraintest is True:
#                myopt.do_roc_train_test()
#            if doplotdistr is True:
#                myopt.do_plot_model_pred()
#            if doimportance is True:
#                myopt.do_importance()
#            if dogridsearch is True:
#                myopt.do_grid()
#            if doboundary is True:
#                myopt.do_boundary()
#            if doefficiencyml is True:
#                myopt.do_efficiency()
#            if dosignifopt is True:
#                myopt.do_significance()
#            if doscancuts is True:
#                myopt.do_scancuts()
#            index = index + 1

    if doml is True and domloption == 2:
        index = 0
        for binmin, binmax in zip(binminarray, binmaxarray):
            myopthipe4ml = Optimiserhipe4ml(data_param[case], case, typean,
                                            binmin, binmax,
                                            training_vars[index],
                                            hipe4ml_hyper_pars[index],
                                            raahp[index])

            if opti_hyperpar_hipe4ml is True:
                myopthipe4ml.do_hipe4mlhyperparopti()
            else:
                myopthipe4ml.set_hipe4ml_modelpar()

            if dotraining:
                myopthipe4ml.do_hipe4mltrain()
                myopthipe4ml.do_hipe4mlplot()
            elif domlplots is True and dotraining is False:
                myopthipe4ml.do_hipe4mlplot(False)
            else:
                myopthipe4ml.get_hipe4mlmodel()
            if dosignifopt is True:
                myopthipe4ml.do_significance(index)
            index = index + 1

    if doml is True and domloption == 3:
        index = 0
        for binmin, binmax in zip([2], [3]):
            myopthipe4ml = Optimiserhipe4mltree(data_param[case], binmin, binmax,
                                                training_vars[index],
                                                hipe4ml_hyper_pars[index])

            if opti_hyperpar_hipe4ml is True:
                myopthipe4ml.do_hipe4mlhyperparopti()
            else:
                myopthipe4ml.set_hipe4ml_modelpar()

            myopthipe4ml.do_hipe4mltrain()
            myopthipe4ml.do_hipe4mlplot()
            index = index + 1

def load_config(user_path: str, default_path: tuple) -> dict:
    """
    Quickly extract either configuration given by user and fall back to package default if no user
    config given.
    Args:
        user_path: path to YAML file
        default_path: tuple were to find the resource and name of resource
    Returns:
        dictionary built from YAML
    """
    logger = get_logger()
    stream = None
    if user_path is None:
        print(default_path[0], default_path[1])
        stream = resource_stream(default_path[0], default_path[1])
    else:
        if not exists(user_path):
            logger_string = f"The file {user_path} does not exist."
            logger.fatal(logger_string)
        stream = open(user_path)
    return yaml.load(stream, yaml.FullLoader)

def main():
    """
    This is used as the entry point for ml-analysis.
    Read optional command line arguments and launch the analysis.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", help="activate debug log level")
    parser.add_argument("--log-file", dest="log_file", help="file to print the log to")
    parser.add_argument("--run-config", "-r", dest="run_config",
                        help="the run configuration to be used")
    parser.add_argument("--database-analysis", "-d", dest="database_analysis",
                        help="analysis database to be used")
    parser.add_argument("--database-ml-models", dest="database_ml_models",
                        help="ml model database to be used")
    parser.add_argument("--database-run-list", dest="database_run_list",
                        help="run list database to be used")
    parser.add_argument("--analysis", "-a", dest="type_ana",
                        help="choose type of analysis")

    args = parser.parse_args()

    configure_logger(args.debug, args.log_file)

    # Extract which database and run config to be used
    pkg_data = "database"
    pkg_data_run_config = "MLTraining"
    run_config = load_config(args.run_config, (pkg_data_run_config, "default_complete.yml"))
    case = run_config["case"]
    if args.type_ana is not None:
        run_config["analysis"]["type"] = args.type_ana

    db_analysis_default_name = f"database_ml_parameters_{case}.yml"
    print(args.database_analysis)
    db_analysis = load_config(args.database_analysis, (pkg_data, db_analysis_default_name))
    db_ml_models = load_config("database/config_model_parameters.yml", (pkg_data, "config_model_parameters.yml"))

    # Run the chain
    do_training(run_config, db_analysis, db_ml_models)
