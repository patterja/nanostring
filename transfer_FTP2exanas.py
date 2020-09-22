#!/usr/bin/env python

## USAGE:
#   transfer_FTP2exanas.py --source_host host --source_user user --source_pass password --dest destination dirname
import argparse
import os
import json
from nanostring_client import nCounter

VERSION="1.0.0"


def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='Transfer to exanas')
    parser.add_argument('--config', type=str, help='nanostring_config.json', default='nanostring_config.json')

    args = parser.parse_args()
    return args



def main():
    args = supply_args()
    with open(args.config, 'r') as cfgfile:
        cfg = json.load(cfgfile)

    #Back up entire nCounter system to exanas01
    ftp_host = nCounter(cfg['nCounterFTP']['host'], cfg['nCounterFTP']['user'], cfg['nCounterFTP']['password'])
    backdirs = ["RCCData", "RLFData", "CDFData"]
    for dir in backdirs:
        destdir = os.path.join(cfg['exanas']['exanas_destdir'], dir)
        if os.path.isdir(destdir):
            ftp_host.download_datadir(os.path.join("/",cfg['nCounterFTP']['user']),destdir)
        else:
            os.mkdir(destdir)
            ftp_host.download_datadir(os.path.join("/", cfg['nCounterFTP']['user'], dir), destdir)
        print("Backed up ", dir)

    ftp_host.close()

if __name__ == "__main__":
    main()
