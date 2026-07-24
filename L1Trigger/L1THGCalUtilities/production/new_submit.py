#!/usr/bin/env python3

import os
import sys
import yaml
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser
from http.client import HTTPException
from multiprocessing import Process
import copy

from CRABClient.UserUtilities import config, ClientException
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException

sys.path.append(".")
from production_tag import production_tag

requestname_base = "eertorer"
output_site = "T3_US_FNALLPC"
output_lfn_base = "/store/user/eertorer/{production_tag}".format(
    production_tag=production_tag
)

def submit(myconfig):
    print("DEBUG : In submit()")
    try:
        crabCommand('submit', config=myconfig)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
        print(hte)
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-y', '--yaml', default='samples_datatest.yaml',
                        help='File with dataset descriptions')
    args = parser.parse_args()

    with open(args.yaml) as f:
        doc = yaml.load(f, Loader=yaml.SafeLoader)
        defaults = doc['defaults'] if 'defaults' in doc else {}

        for sample in sorted(doc["samples"].keys()):
            info = copy.deepcopy(defaults)
            info.update(doc["samples"][sample])
            print("\n\n*** Sample {} ***".format(sample))

            for dataset_shortname, dataset in info['datasets'].items():
                print("\n*** Submitting {}: {}".format(dataset_shortname, dataset))

                isMC = info.get("isMC", None)
                if isMC is None:
                    raise ValueError("Please specify parameter isMC in the YAML.")

                this_config = config()

                # ---- General ----
                this_config.section_('General')
                this_config.General.transferOutputs = True
                this_config.General.transferLogs = True
                this_config.General.workArea = "crab/{}_{}/".format(requestname_base, production_tag)
                this_config.General.requestName = "{}_{}_{}_{}".format(
                    requestname_base, production_tag, info["year"], dataset_shortname
                )

                # ---- JobType ----
                this_config.section_('JobType')
                this_config.JobType.pluginName = 'Analysis'
                this_config.JobType.psetName = os.path.expandvars(info.get("pset", None))
                this_config.JobType.allowUndistributedCMSSW = True
                this_config.JobType.numCores = 4
                this_config.JobType.maxMemoryMB = 8000
                globaltag = info.get("globaltag", None)
                this_config.JobType.pyCfgParams = [
                    f'isMC={isMC}',
                    'reportEvery=1000',
                    f'tag={production_tag}',
                    f'globalTag={globaltag}'
                ]

                # ---- Site ----
                this_config.section_('Site')
                this_config.Site.storageSite = output_site

                # ---- Data ----
                this_config.section_('Data')
                this_config.Data.publication = False
                this_config.Data.outLFNDirBase = f"{output_lfn_base}/{info['year']}/{sample}"
                this_config.Data.outputDatasetTag = dataset_shortname
                this_config.Data.inputDBS = 'global'

                # Handle either a text file (userInputFiles) or a DBS dataset
                if dataset.endswith(".txt"):
                    # We assume this is a file list
                    print("INFO: Detected a file list (not a DBS dataset). Using userInputFiles.")
                    file_list = []
                    with open(dataset) as fileobj:
                        for line in fileobj:
                            line = line.strip()
                            if line:
                                file_list.append(line)
                    
                    this_config.Data.userInputFiles = file_list
                    # Force file-based splitting
                    this_config.Data.splitting = "FileBased"

                    # If 'unitsPerJob' is in the YAML, use it. Otherwise default to 1.
                    unitsPerJob = info.get("unitsPerJob", 1)
                    this_config.Data.unitsPerJob = unitsPerJob

                else:
                    # We assume this is a DBS dataset
                    print("INFO: Detected a normal DBS dataset. Using inputDataset.")
                    this_config.Data.inputDataset = dataset

                    splitting_mode = info.get("splitting", "Automatic")
                    if splitting_mode not in ["Automatic", "FileBased", "LumiBased"]:
                        raise ValueError(f"Unrecognized splitting mode: {splitting_mode}")
                    this_config.Data.splitting = splitting_mode

                    # If unitsPerJob/totalUnits are set, apply them
                    unitsPerJob = info.get("unitsPerJob", None)
                    if unitsPerJob is not None:
                        this_config.Data.unitsPerJob = unitsPerJob

                    totalUnits = info.get("totalUnits", None)
                    if totalUnits is not None:
                        this_config.Data.totalUnits = totalUnits

                # If real data, load lumimask
                if not isMC:
                    this_config.Data.lumiMask = info.get('lumimask', None)
                else:
                    this_config.Data.lumiMask = ''

                allowInvalid = info.get("allowInvalid", False)
                if allowInvalid:
                    this_config.Data.allowNonValidInputDataset = True

                print(this_config)
                p = Process(target=submit, args=(this_config,))
                p.start()
                p.join()

            print("*** Done with Sample {} ***\n\n".format(sample))
