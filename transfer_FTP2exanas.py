#!/usr/bin/env python

## USAGE:
#   transfer_FTP2exanas.py --source_host host --source_user user --source_pass password --dest destination dirname
import argparse
import os
from nanostring_client import nCounter

VERSION="1.0.0"


def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='Transfer batch from ftp to mounted X drive \n'
                                                 'Remote = remote directory \n'
                                                 'Source = FTP Nanostring machine')
    parser.add_argument('--source_host', type=str, help='host name')
    parser.add_argument('--source_user', type=str, help='user name')
    parser.add_argument('--source_pass', type=str, help='password')
    parser.add_argument('--dest', type=str, help='remote automated_data dir path',
                        default="/home/exacloud/clinicalarchive/KDL/Nanostring")
    args = parser.parse_args()
    return args



def main():
    args = supply_args()
    ftp_host = nCounter(args.source_host, args.source_user, args.source_pass)
    backdirs = ["RCCData", "RLFData", "CDFData"]
    for dir in backdirs:
        destdir = os.path.join(args.dest, dir)
        if os.path.isdir(destdir):
            ftp_host.download_datadir(os.path.join("/",args.source_user, dir),destdir)
        else:
            os.mkdir(destdir)
            ftp_host.download_datadir(os.path.join("/", args.source_user, dir), destdir)
        print("Backed up ", dir)

    ftp_host.close()


if __name__ == "__main__":
    main()
