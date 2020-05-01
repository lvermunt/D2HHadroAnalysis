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

from CreateProcessTTree.multiprocesser import MultiProcesser
from CreateProcessTTree.processer import Processer
from machine_learning_hep.logger import configure_logger, get_logger
from machine_learning_hep.utilities import checkdirlist
from machine_learning_hep.utilities import checkmakedirlist

try:
    import logging
    import absl.logging
    logging.root.removeHandler(absl.logging._absl_handler) # pylint: disable=protected-access
    absl.logging._warn_preinit_stderr = False # pylint: disable=protected-access
except Exception as e: # pylint: disable=broad-except
    print("##############################")
    print("Failed to fix absl logging bug", e)
    print("##############################")


def do_application(data_config: dict, data_param: dict, run_param: dict):

    # Disable any graphical stuff. No TCanvases opened and shown by default
    gROOT.SetBatch(True)

    logger = get_logger()
    logger.info("Do analysis chain")

    # If we are here we are interested in the very first key in the parameters database
    for k in data_param.keys():
        case = k
        break

    doapplydata = data_config["mlapplication"]["data"]["doapply"]
    doapplymc = data_config["mlapplication"]["mc"]["doapply"]
    doapplydatamodelhandler = data_config["mlapplication"]["data"]["doapplymodelhandler"]
    doapplymcmodelhandler = data_config["mlapplication"]["mc"]["doapplymodelhandler"]
    domergeapplydata = data_config["mlapplication"]["data"]["domergeapply"]
    domergeapplymc = data_config["mlapplication"]["mc"]["domergeapply"]
    docontinueapplydata = data_config["mlapplication"]["data"]["docontinueafterstop"]
    docontinueapplymc = data_config["mlapplication"]["mc"]["docontinueafterstop"]

    domlprefilterstep = data_config["mlapplication"]["domlprefilterstep"]
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

    dirpklskdecmc = data_param[case]["mlapplication"]["mc"]["pkl_skimmed_dec"]
    dirpklskdec_mergedmc = data_param[case]["mlapplication"]["mc"]["pkl_skimmed_decmerged"]
    dirpklskdecdata = data_param[case]["mlapplication"]["data"]["pkl_skimmed_dec"]
    dirpklskdec_mergeddata = data_param[case]["mlapplication"]["data"]["pkl_skimmed_decmerged"]
    if domlprefilterstep_indb is True:
        dirpklskdecmc = [d + "/prefilter" for d in dirpklskdecmc]
        dirpklskdec_mergedmc = [d + "/prefilter" for d in dirpklskdec_mergedmc]
        dirpklskdecdata = [d + "/prefilter" for d in dirpklskdecdata]
        dirpklskdec_mergeddata = [d + "/prefilter" for d in dirpklskdec_mergeddata]
        data_param[case]["mlapplication"]["mc"]["pkl_skimmed_dec"] = dirpklskdecmc
        data_param[case]["mlapplication"]["mc"]["pkl_skimmed_decmerged"] = dirpklskdec_mergedmc
        data_param[case]["mlapplication"]["data"]["pkl_skimmed_dec"] = dirpklskdecdata
        data_param[case]["mlapplication"]["data"]["pkl_skimmed_decmerged"] = dirpklskdec_mergeddata
    if domlprefilterstep_indb is False:
        dirpklskdecmc = [d + "/analysis" for d in dirpklskdecmc]
        dirpklskdec_mergedmc = [d + "/analysis" for d in dirpklskdec_mergedmc]
        dirpklskdecdata = [d + "/analysis" for d in dirpklskdecdata]
        dirpklskdec_mergeddata = [d + "/analysis" for d in dirpklskdec_mergeddata]
        data_param[case]["mlapplication"]["mc"]["pkl_skimmed_dec"] = dirpklskdecmc
        data_param[case]["mlapplication"]["mc"]["pkl_skimmed_decmerged"] = dirpklskdec_mergedmc
        data_param[case]["mlapplication"]["data"]["pkl_skimmed_dec"] = dirpklskdecdata
        data_param[case]["mlapplication"]["data"]["pkl_skimmed_decmerged"] = dirpklskdec_mergeddata

    #creating folder if not present
    counter = 0
    if docontinueapplymc is False:
        if doapplymc is True:
            counter = counter + checkdirlist(dirpklskdecmc)

        if domergeapplymc is True:
            counter = counter + checkdirlist(dirpklskdec_mergedmc)

    if docontinueapplydata is False:
        if doapplydata is True:
            counter = counter + checkdirlist(dirpklskdecdata)

        if domergeapplydata is True:
            counter = counter + checkdirlist(dirpklskdec_mergeddata)

    if counter < 0:
        sys.exit()

    # check and create directories
    if docontinueapplymc is False:
        if doapplymc is True:
            checkmakedirlist(dirpklskdecmc)

        if domergeapplymc is True:
            checkmakedirlist(dirpklskdec_mergedmc)

    if docontinueapplydata is False:
        if doapplydata is True:
            checkmakedirlist(dirpklskdecdata)

        if domergeapplydata is True:
            checkmakedirlist(dirpklskdec_mergeddata)

    proc_class = Processer
    mymultiprocessmc = MultiProcesser(case, proc_class, data_param[case], run_param, "mc")
    mymultiprocessdata = MultiProcesser(case, proc_class, data_param[case], run_param, "data")

    #perform the analysis flow
    if doapplydata is True:
        if doapplydatamodelhandler is True:
            mymultiprocessdata.multi_apply_hipe4ml_allperiods()
        else:
            mymultiprocessdata.multi_apply_allperiods()
    if doapplymc is True:
        if doapplydatamodelhandler is True:
            mymultiprocessmc.multi_apply_hipe4ml_allperiods()
        else:
            mymultiprocessmc.multi_apply_allperiods()
    if domergeapplydata is True:
        mymultiprocessdata.multi_mergeapply_allperiods()
    if domergeapplymc is True:
        mymultiprocessmc.multi_mergeapply_allperiods()

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
    db_run_list = load_config("database/database_run_list.yml", (pkg_data, "database_run_list.yml"))

    # Run the chain
    do_application(run_config, db_analysis, db_run_list)
