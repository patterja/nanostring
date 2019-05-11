#!/usr/bin/env python

#USAGE:
# Nanostring output RCC readouts (html format) to raw data matrix
#NOTES:
#TODO: break out test: checking sample metrics and check for unique sample names in sample sheet
#ARGUMENTS:
#RETURN:
import os
import re
import argparse
import pandas
import xml.etree.ElementTree as ET

VERSION="0.1.0"



def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='Nanostring output RCC readouts (html format) and '
                                                 'converts to matrix, geomean normalized values')
    parser.add_argument('samplesheet', help='samplesheet.txt')
    parser.add_argument('rcc_dir', help='raw RCC files directory')
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
        df_samp_attrib = pandas.DataFrame(header, samp_attrib).T
        df_lane_attrib = pandas.DataFrame(lane_attrib[1:len(lane_attrib)], columns=lane_attrib[0])

    return(dfcounts, df_samp_attrib, df_lane_attrib)


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
    with open(samplesheet, 'r') as ss_handle:
        for line in ss_handle:
            items = line.strip().split("\t")
            ss_dict[items[0]] = items[1]
    return(ss_dict)


def main():

    args = supply_args()

    sampleid = parse_samplesheet(args.samplesheet)
    rcc_counts_dict = {}
    samp_attrib_dict = {}
    lane_attrib_dict = {}

    for file in os.listdir(args.rcc_dir):
        rcc_dir = args.rcc_dir
        if file.endswith(".RCC"):

            #get sample number
            samp_number = re.sub(".RCC", "", file.split("_")[-1])
            dfrcc, df_samp_attrib, df_lane_attrib = parseRCC(os.path.join(rcc_dir + "/", file))

            #rename column name with sample id
            dfrcc.rename(columns={'Count': sampleid[file]}, inplace=True)
            rcc_counts_dict[sampleid[file]] = dfrcc

            df_lane_attrib.rename(columns={df_lane_attrib.columns[1]: sampleid[file]}, inplace=True)
            lane_attrib_dict[df_lane_attrib.columns[1]] = df_lane_attrib

            samp_attrib_dict[sampleid[file]] = df_samp_attrib

            batch_name=re.sub(r'_%s+.RCC' % samp_number, "", file)

    #merge dictionary of dataframes together
    raw_data = reduce(lambda x, y: pandas.merge(x, y, on=['CodeClass', 'Name', 'Accession']),
                      rcc_counts_dict.values())

    lane_attrib = reduce(lambda x, y: pandas.merge(x, y, on=['ID']),
                         lane_attrib_dict.values())

    raw_data.to_csv(args.rcc_dir + "/"+ batch_name+ "_rawdata.txt", sep='\t', index=False)

    #test to see if samp_attribs are all the same
    for v in range(len(samp_attrib_dict.values()) - 1):
        if samp_attrib_dict.values()[v].equals(samp_attrib_dict.values()[v+1]) != True:
            print("Samples are not from one batch. Sample Attributes differ")
            print(samp_attrib_dict.keys()[v])
        else:
            print("Samples are from one batch. OK")

    with open(args.rcc_dir + "/" + batch_name + "_run_metrics.txt", 'w') as met:
        met.write(samp_attrib_dict.values()[0].to_csv(sep="\t"))
        met.write(lane_attrib.to_csv(sep="\t"))


if __name__ == "__main__":
    main()
