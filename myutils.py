import os
import sys
from typing import List, Set, Dict, Tuple

STUB = False

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class Logger:
    def __init__(self, logname: str = None) -> None:
        self.logname: str = logname
    def set_file(self, logname: str) -> None:
        self.logname = logname
    def log2file(self, msg: str):
        if self.logname is not None:
            with open(self.logname, "a") as fo:
                fo.write(str(msg) + '\n')
    def log(self, msg: str, col: bcolors = bcolors.BOLD) -> None:
        print(col + str(msg) + bcolors.ENDC)
        self.log2file(msg)

logger = Logger()

def create_directories(dirlist: List[str]) -> None:
    logger.log("#Creating directories", bcolors.BOLD)
    for d in dirlist:
        if not os.path.isdir(d):
            try:
                os.mkdir(d)
                logger.log(f"#----{d} created", bcolors.OKGREEN)
            except Exception as e:
                logger.log(e, bcolors.FAIL)
                logger.log(f"#----Unable to create directory: {d}", bcolors.FAIL)
        else:
            logger.log(f"#----{d} already exists", bcolors.WARNING)

def run_command(cmd: str):
    logger.log(f"#----Running command:\n{cmd}", bcolors.HEADER)
    if STUB:
        logger.log(f"#---STUB mode ON: Skipping", bcolors.OKGREEN)
        return
    try:
        exitcode: int = os.system(cmd)
        if exitcode == 0:
            logger.log(f"#---Success: {cmd}", bcolors.OKGREEN)
        else:
            logger.log(f"#---FAIL {exitcode}: {cmd}", bcolors.FAIL)
    except Exception as e:
        logger.log(e, bcolors.FAIL)
        logger.log(f"#---FAIL: {cmd}", bcolors.FAIL)