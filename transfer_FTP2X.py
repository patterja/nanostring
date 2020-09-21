#!/usr/bin/env python

## USAGE:
#   transfer_nanostring.py --source_host --source_user --source_pass --batch --remote_data_dir --remote_ss_dir

## ARGUMENTS:
#
## RETURN:
#
import argparse
import os
import sys
import shutil
import zipfile
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
    parser.add_argument('--batch', type=str, help='batch id')
    parser.add_argument('--dest_data_dir', type=str, help='dest automated_data dir path')
    parser.add_argument('--dest_ss_dir', type=str, help='dest samplesheet dir path')
    args = parser.parse_args()
    return args

def prep_dest_dir(batch, dest_data_dir, dest_ss_dir):
    ss = os.path.join(dest_ss_dir, batch +"_samplesheet.txt")
    data_dir = os.path.join(dest_data_dir, batch)

    if os.path.isdir(data_dir):
        print("Directory exists")
    else:
        os.mkdir(os.path.join(dest_data_dir, batch))

    if os.path.exists(ss):
        shutil.copy(ss, data_dir, follow_symlinks=True)
    else:
        print("No samplesheet exists for " + batch)




def main():
    args = supply_args()
    zfile = args.batch + "_RCC.ZIP"
    data_dir = os.path.join(args.dest_data_dir, args.batch)
    zfile_path = os.path.join(data_dir, zfile)
    ftp_host = nCounter(args.source_host, args.source_user, args.source_pass)
    #ftp_host = ftputil.FTPHost("10.125.46.25", "technician", "MAX123")
    # check if zip file and directory structure already exist in automated data
    if os.path.exists(data_dir):
        print(os.path.join(args.batch + "already exists"))
    else:
        # create automated_data directory and copy samplesheet
        prep_dest_dir(args.batch, args.dest_data_dir, args.dest_ss_dir)
        ftp_host.download_batch(args.batch, data_dir)

    if os.path.isfile(zfile_path):
        with zipfile.ZipFile(zfile_path,'r') as zobj:
            zobj.extractall(data_dir)
    #move zip file to zipfile directory
    # shutil.move(zfile_path, os.path.join(args.dest_data_dir, "RCCData", zfile))
    os.remove(zfile_path)


if __name__ == "__main__":
    main()
