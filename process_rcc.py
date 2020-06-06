#!/usr/bin/env python

## USAGE:
#   Nanostring output RCC readouts (html format) to raw data matrix
## NOTES:
#   TODO: break out tests
#       checking sample metrics
#       check for unique sample names in sample sheet
#
## ARGUMENTS:
#   samplesheet: tsv of RCC filename in one col and sample name, 'ANTIBODY_REFERENCE.csv'
#   rcc files
#   abfile: antibody file to rename antibody names.
## RETURN:
#   rawdata.txt: tab-sep file with table of Antibody x Sample
#   run_metrics.txt
#
import re
import argparse
import pandas
import xml.etree.ElementTree as ET
from functools import reduce

VERSION = "2.0.0"


def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='Nanostring output RCC readouts (html format) and '
                                                 'converts to raw data tsv')
    parser.add_argument('rcc_files', type=str, nargs='+', help='raw RCC files')
    parser.add_argument('--samplesheet', type=str, help='samplesheet.txt')
    parser.add_argument('--abfile', type=argparse.FileType('r'), help='ANTIBODY_REFERENCE.csv')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def parseRCC(rcc_file):
    """
    Take RCC file output df
    Args:
        path to RCC
    Returns:
        dataframe
    Examples:
    """
    dict_rcc = {}
    with open(rcc_file, 'r') as xml_handle:
        complete_xml = "<root>" + str(xml_handle.read()) + "</root>"
        root = ET.fromstring(complete_xml)
        for child in root:
            dict_rcc[child.tag] = child.text

        counts = [line.strip().split(',') for line in dict_rcc['Code_Summary'].strip().split('\n')]
        header = [line.strip().split(',') for line in dict_rcc['Header'].strip().split('\n')]
        samp_attrib = [line.strip().split(',') for line in
                       dict_rcc['Sample_Attributes'].strip().split('\n')]
        lane_attrib = [line.strip().split(',') for line in
                       dict_rcc['Lane_Attributes'].strip().split('\n')]

        dfcounts = pandas.DataFrame(counts[1:len(counts)], columns=counts[0])
        df_lane_attrib = pandas.DataFrame(lane_attrib[1:len(lane_attrib)], columns=lane_attrib[0])
        # df_samp_attrib = pandas.DataFrame(header, samp_attrib).T
        # df_samp_attrib = pandas.DataFrame(header[0:1], samp_attrib).T
    return (dfcounts, header, samp_attrib, df_lane_attrib)

def parse_Ab_ref(abfile):
    """
    get antibody shortnames.
    Args:
        ab_handle
    Returns:
        dict
    Examples:
    """
    ab_name ={}
    mouseAb = []
    rabbitAb =[]
    for line in abfile:
        if not line.lstrip().startswith('#'):
            items = line.strip().split(",")
            ab_name[items[1]] = items[0]
            if items[4] == "mouse":
                mouseAb.append(items[0])
            else:
                rabbitAb.append(items[0])
    return(ab_name, rabbitAb, mouseAb)


def parse_samplesheet(samplesheet):
    """
    get sample names, expecting to get more complex with this.
    ie. make directory structure to store data, etc
    Args:
        samplesheet
    Returns:
        dict
    Examples:
    """
    ss_dict = {}
    with open(samplesheet, 'r') as ss:
        count = 0
        for line in ss:
            items = line.strip().split("\t")
            ss_dict[count] = items[1]
            count += 1
    return (ss_dict)


def main():
    args = supply_args()
    if len(args.rcc_files) < 12:
        lrcc = len(args.rcc_files)
        print("\nOnly " + str(lrcc) + " RCC files specified!\nContinuing with " + str(
            lrcc) + " files.\n")

    sampleid = parse_samplesheet(args.samplesheet)
    print(sampleid)

    rcc_counts_dict = {}
    samp_attrib_dict = {}
    header_dict = {}
    lane_attrib_dict = {}


    for file in args.rcc_files:
        if file.endswith(".RCC"):
            # get sample number
            # samp_number = re.sub(".RCC", "", file.split("_")[-1])
            dfrcc, header, samp_attrib, df_lane_attrib = parseRCC(file)

            samp_number = int(df_lane_attrib.columns[1])
            # rename column name with sample id for rcc files
            dfrcc.rename(columns={'Count': sampleid[samp_number]}, inplace=True)
            rcc_counts_dict[sampleid[samp_number]] = dfrcc

            # df_lane_attrib.rename(columns={df_lane_attrib.columns[1]: sampleid[file]}, inplace=True)
            df_lane_attrib = df_lane_attrib.append(
                pandas.Series(['SampleName', sampleid[samp_number]],
                              index=df_lane_attrib.columns), ignore_index=True)
            lane_attrib_dict[df_lane_attrib.columns[1]] = df_lane_attrib

            samp_attrib_dict[sampleid[samp_number]] = samp_attrib
            header_dict[sampleid[samp_number]] = header

    # Check header and samp_attribs

    for idx in range(1, len(header_dict)):
        if list(header_dict.values())[idx - 1] != list(header_dict.values())[idx]:
            print("RCC header are not equal")
            # print(header_dict.values()[idx])

    # test to see if samp_attribs are all the same, all RCC files are from same run
    for idx in range(1, len(samp_attrib_dict)):
        if list(samp_attrib_dict.values())[idx - 1] != list(samp_attrib_dict.values())[idx]:
            print("RCC Sample Attributes are not equal")
            # print(samp_attrib_dict.values()[idx])
            if samp_attrib_dict.values()[idx][0][1:] == samp_attrib_dict.values()[idx][0][1:]:
                print("Actually only ID don't match, is this an old RCC file or have chaned RLFs")
            else:
                print("RCC Sample Attributes IDs and values are not equal. Stopping")
                print("Samples are not from one batch. Sample Attributes differ")
                exit("Error, Check your RCC files")
        else:
            print("Samples attribs are okay and are from one batch. OK")

    # merge dictionary of dataframes together
    raw_data = reduce(lambda x, y: pandas.merge(x, y, on=['CodeClass', 'Name', 'Accession']),
                      rcc_counts_dict.values())

    lane_attrib_combined = reduce(lambda x, y: pandas.merge(x, y, on=['ID']),
                                  lane_attrib_dict.values())

    with open("run_metrics.txt", 'w') as met:
        for item in list(header_dict.values())[0]:
            met.write('\t'.join(item))
            met.write('\n')
        for item in list(samp_attrib_dict.values())[0]:
            met.write('\t'.join(item))
            met.write('\n')
        met.write(lane_attrib_combined.to_csv(sep="\t"))

    # Change long name to something else if necessary
    with args.abfile:
        ab_name, rabbitAb, mouseAb = parse_Ab_ref(args.abfile)

    # Write raw data file
    for name in raw_data["Name"].iteritems():
        if not re.search("POS|NEG", name[1]):
            # print(name)
            raw_data['Name'][name[0]] = ab_name[name[1].strip().split("|")[0]]

    raw_data.to_csv("rawdata.txt", sep='\t', index=False)

if __name__ == "__main__":
    main()
