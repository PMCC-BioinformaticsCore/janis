import argparse

import dummywf;
import json;
import os;
import uuid;
import logging;

def main( opts ):
    print("WeHi Pipeline Definition Demo!");

    sessionID = str(uuid.uuid1());
    print("Session: " + sessionID);

    configFile = opts.config;
    if configFile is None:
        configFile = "config.json";
    else:
        configFile = str(configFile);

    configPath = os.path.abspath(configFile);
    print("Using config file: " + configPath);

    conf = None;
    if os.path.isfile(configPath):
        confFile = open(configPath);
        conf = json.load(confFile);
    else:
        print("Config file does not exists.");

    try:
        wsDir = conf["workspace"];
    except:
        wsDir = None;

    if wsDir is None:
        wsDir = ".";
    else:
        wsDir = str(wsDir);
    wsPath = os.path.abspath(wsDir);
    print("Using Workspace: " + wsPath );

    try:
        logPrefix = conf["logfileprefix"];
    except:
        logPrefix = None;

    if logPrefix is None:
        logPrefix = "wehi-util.log";

    logFile = wsPath + "/" + logPrefix + "-" + sessionID + ".txt";
    logPath = os.path.abspath(logFile);

    print("Using log file: " + logPath );

    logger = logging.getLogger("WEHIPLUTIL");
    logHandler = logging.FileHandler(logPath);
    formatter = logging.Formatter("%(message)s");
    logHandler.setFormatter(formatter);
    logger.addHandler( logHandler );
    logger.setLevel( logging.INFO );

    cwlDepot = None;
    try:
        cwlDepot = conf["cwl.depot"];
    except:
        cwlDepot = None;

    if cwlDepot is None:
        raise AssertionError("CWL Depot path not configured.");
    else:
        cwlDepotPath = str(cwlDepot);
    cwlDepotPath = os.path.abspath(cwlDepotPath);
    logger.info("Using CWL Depot: " + cwlDepotPath );

    cwlPath = wsPath + "/dummy.cwl";
    cwlPath = os.path.abspath(cwlPath);
    logger.info("Using CWL file: " + str(cwlPath));

    wf = dummywf.DummyWorkflow();
    wf.emitCWL( cwlPath, logger );

if __name__ == "__main__":
    argprsr = argparse.ArgumentParser();
    argprsr.add_argument("--config", required=False);
    opts = argprsr.parse_args()
    main( opts );


#Python exceptions
#https://docs.python.org/2/library/exceptions.html#exception-hierarchy

