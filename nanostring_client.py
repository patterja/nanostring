#!/usr/bin/env python

## author: Janice Patterson
## targets nCounter MAX/FLEX system

#import ftputil
from ftputil import FTPHost
import sys
import os

class nCounter(FTPHost):
    '''
    inheriting from ftputil.FTPHost
    '''

    def download(self, source_dir, dest_dir):
        '''
        :param source_dir: /technician/RCCData, /technician/RLFData, /technician/CLFData
        :param dest_dir: dest directory name the directory RCCData, RLFData, CLFData
        :return: transfers
        '''
        if self.path.exists(source_dir):
            recursive = self.walk(source_dir, topdown=True, onerror=None)
            for root, dirs, files in recursive:
                for file in files:
                    print(root + "/" + file + " to ", dest_dir)
                    fpath = self.path.join(root, file)
                    print(fpath)
                    if self.path.isfile(fpath):
                        dest_file = os.path.join(dest_dir, file)
                        print(dest_file)
                        # download only if newer double make sure
                        if not self.path.exists(dest_file):
                            print("does not exist")
                            self.download_if_newer(fpath,dest_file)
                        else:
                            print(dest_file + " exists already, skipping.")
        else:
            print("FTP dir is "+ self.getcwd())
            print(self.listdir(self.curdir))

    def download_batch(self, batch, user):
        '''
        :param batch: batch id only, path to file /technician/RCCData is taken care of
        :return: download zip to current directory
        '''
        zfile = batch + "_RCC.ZIP"
        batchpath = self.path.join("/", user, "RCCData")
        zfiles = self.listdir(batchpath)
        if zfile in zfiles:
            self.download(os.path.join(batchpath, zfile), zfile)
        else:
            print("no " + zfile)
            print(zfiles)


