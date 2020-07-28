#!/usr/bin/env python

## USAGE:
#   transfer_nanostring.py --source_host --source_user --source_pass --batch --remote_data_dir --remote_ss_dir

## ARGUMENTS:
#
## RETURN:
#
import argparse
import ftputil
import os
import sys
import shutil
import zipfile

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
    parser.add_argument('--remote_data_dir', type=str, help='remote automated_data dir path')
    parser.add_argument('--remote_ss_dir', type=str, help='remote samplesheet dir path')
    args = parser.parse_args()
    return args

def prep_remote_dir(batch, remote_data_dir, remote_ss_dir):
    ss = os.path.join(remote_ss_dir, batch +"_samplesheet.txt")
    data_dir = os.path.join(remote_data_dir, batch)

    if os.path.isdir(data_dir):
        print("Directory exists")
    else:
        os.mkdir(os.path.join(remote_data_dir, batch))

    if os.path.exists(ss):
        shutil.copy(ss, data_dir, follow_symlinks=True)
    else:
        sys.exit("")



def download(host, user, password):
    # Download some files from the login directory.
    with ftputil.FTPHost(host, user, password) as ftp_host:
        names = ftp_host.listdir(ftp_host.curdir)
        for name in names:
            if ftp_host.path.isfile(name):
                print(name)
                ftp_host.download(name, name)  # remote, local
        # Make a new directory and copy a remote file into it.
        ftp_host.mkdir("newdir")
        with ftp_host.open("index.html", "rb") as source:
            with ftp_host.open("newdir/index.html", "wb") as target:
                ftp_host.copyfileobj(source, target)  # similar to shutil.copyfileobj
def main():
    args = supply_args()
    zfile = args.batch + "_RCC.ZIP"
    zfile_path = os.path.join(args.remote_data_dir, zfile)
    data_dir = os.path.join(args.remote_data_dir, args.batch)

    #ftp_host = ftputil.FTPHost("10.125.46.25", "technician", "MAX123")
    # check if zip file and directory structure already exist in automated data
    if os.path.exists(zfile_path):
        print(os.path.join(zfile_path + " exists"))
    else:
        with ftputil.FTPHost(args.source_host, args.source_user, args.source_pass) as ftp_host:
            if zfile in ftp_host.listdir("./RCCData"):
                ftp_host.download(os.path.join("./RCCData", zfile), zfile_path)

    #create automated_data directory and copy samplesheet
    prep_remote_dir(args.batch, args.remote_data_dir, args.remote_ss_dir)

    if os.path.isdir(data_dir):
        with zipfile.ZipFile(zfile_path,'r') as zobj:
            zobj.extractall(data_dir)

if __name__ == "__main__":
    main()
