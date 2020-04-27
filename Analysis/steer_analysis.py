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

import yaml
from pkg_resources import resource_stream
# To set batch mode immediately
from ROOT import gROOT  # pylint: disable=import-error, no-name-in-module

from Analysis.analyser_manager import AnalyserManager
from Analysis.analyseritsupgrade import AnalyserITSUpgrade
from machine_learning_hep.logger import configure_logger, get_logger
from machine_learning_hep.utilities import checkdirlist, checkdir
from machine_learning_hep.utilities import checkmakedirlist, checkmakedir

try:
    import logging
    import absl.logging
    logging.root.removeHandler(absl.logging._absl_handler) # pylint: disable=protected-access
    absl.logging._warn_preinit_stderr = False # pylint: disable=protected-access
except Exception as e: # pylint: disable=broad-except
    print("##############################")
    print("Failed to fix absl logging bug", e)
    print("##############################")


def do_analysis(data_config: dict, data_param: dict): # pylint: disable=too-many-locals, too-many-statements, too-many-branches

    # Disable any graphical stuff. No TCanvases opened and shown by default
    gROOT.SetBatch(True)

    logger = get_logger()
    logger.info("Do analysis chain")

    # If we are here we are interested in the very first key in the parameters database
    for k in data_param.keys():
        case = k
        break

    typean = data_config["analysis"]["type"]
    doanalysismc = data_config["analysis"]["mc"]["activate"]
    doanalysisdata = data_config["analysis"]["data"]["activate"]
    dohistomass = data_config["analysis"]["masshisto"]
    doefficiency = data_config["analysis"]["effhisto"]
    doprobscan = data_config["analysis"]["probscan"]
    dofitbackpars = data_config["analysis"]["fitbackpars"]
    doexpectedsignf = data_config["analysis"]["expectedsignf"]
    dobkgshapestudy = data_config["analysis"]["bkgshapestudy"]

    dirresultsdata = data_param[case]["analysis"][typean]["data"]["results"]
    dirresultsmc = data_param[case]["analysis"][typean]["mc"]["results"]
    dirresultsdatatot = data_param[case]["analysis"][typean]["data"]["resultsallp"]
    dirresultsmctot = data_param[case]["analysis"][typean]["mc"]["resultsallp"]

    #creating folder if not present
    counter = 0

    if doanalysismc is True:
        counter = counter + checkdirlist(dirresultsmc)
        counter = counter + checkdir(dirresultsmctot)

    if doanalysisdata is True:
        counter = counter + checkdirlist(dirresultsdata)
        counter = counter + checkdir(dirresultsdatatot)

    if counter < 0:
        if doprobscan is True:
            sys.exit()
        else:
            logger.warning("Directories already exists (see above), but no new prob scan")
    else:
        # check and create directories
        if doanalysismc is True:
            checkmakedirlist(dirresultsmc)
            checkmakedir(dirresultsmctot)

        if doanalysisdata is True:
            checkmakedirlist(dirresultsdata)
            checkmakedir(dirresultsdatatot)

    ana_class = AnalyserITSUpgrade
    ana_mgr = AnalyserManager(ana_class, data_param[case], case, typean)

    # Collect all desired analysis steps
    analyse_steps = []
    if dohistomass is True:
        analyse_steps.append("invmass")
    if doefficiency is True:
        analyse_steps.append("efficiency")
    if doprobscan is True:
        analyse_steps.append("probability_scan")
    if dofitbackpars is True:
        analyse_steps.append("parametrise_background_scan")
    if doexpectedsignf is True:
        analyse_steps.append("expected_significance_print")
    if dobkgshapestudy is True:
        analyse_steps.append("bkgshapestudy_mass_histo")

    # Now do the analysis
    ana_mgr.analyse(*analyse_steps)

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
    pkg_data_run_config = "CreateProcessTTree"
    run_config = load_config(args.run_config, (pkg_data_run_config, "default_complete.yml"))
    case = run_config["case"]
    if args.type_ana is not None:
        run_config["analysis"]["type"] = args.type_ana

    db_analysis_default_name = f"database_ml_parameters_{case}.yml"
    print(args.database_analysis)
    db_analysis = load_config(args.database_analysis, (pkg_data, db_analysis_default_name))

    # Run the chain
    do_analysis(run_config, db_analysis)
