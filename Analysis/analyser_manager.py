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

from machine_learning_hep.logger import get_logger

# pylint: disable=too-many-instance-attributes
class AnalyserManager:
    """
    Manager class handling analysis and systematic objects
    """

    def __init__(self, ana_class, database, case, typean, *args):

        self.ana_class = ana_class
        self.database = database
        self.case = case
        self.typean = typean

        # Additional arguments to be forwarded to the analysers
        self.add_args = args

        self.logger = get_logger()

        self.analysers = []
        self.after_burner = None

        self.is_initialized = False


    def initialize(self):
        """
        Collect all analyser objects required in a list and initialises the after_burner if present
        """

        if self.is_initialized:
            return

        self.logger.info("Initialize analyser manager for analyser %s", self.ana_class.__name__)

        useperiod = self.database["analysis"][self.typean]["useperiod"]

        for ip, period in enumerate(useperiod):
            self.analysers.append(self.ana_class(self.database, self.case, self.typean, ip,
                                                 *self.add_args))
        self.analysers.append(self.ana_class(self.database, self.case, self.typean, None,
                                             *self.add_args))

        # get after-burner, if any
        self.after_burner = self.analysers[-1].get_after_burner()
        if self.after_burner:
            self.after_burner.analysers = self.analysers[:-1]
            self.after_burner.analyser_merged = self.analysers[-1]

        self.is_initialized = True


    def analyse(self, *ana_steps):
        """
        Gives a list of analysers and analysis steps do each step for each analyser
        Args:
            ana_steps: list of analysis steps as strings
        """

        if not ana_steps:
            self.logger.info("No analysis steps to be done for Analyser class %s. Return...",
                             self.ana_class.__name__)
            return

        self.initialize()

        self.logger.info("Run all registered analysers of type %s for following analysis steps",
                         self.ana_class.__name__)
        for step in ana_steps:
            print(f"  -> {step}")

        # Collect potentially failed systematic steps
        failed_steps = []
        failed_steps_after_burner = []
        for step in ana_steps:
            for analyser in self.analysers[:-1]:
                if not analyser.step(step):
                    failed_steps.append((analyser.__class__.__name__, step))
                    # If analysis step could not be found here,
                    # we don't need to go on trying this steps since all analysers are of the
                    # same class
                    break

            # Run after-burner if one was provided by the analyser object
            if self.after_burner and not self.after_burner.step(step):
                failed_steps_after_burner.append((self.after_burner.__class__.__name__, step))

            # Do analysis step for period-merged analyser
            self.analysers[-1].step(step)

        if failed_steps:
            self.logger.error("Following analysis steps could not be found:")
            for fs in failed_steps:
                print(f"Analyser class: {fs[0]}, analysis step: {fs[1]}")
