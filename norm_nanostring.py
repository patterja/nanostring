#!/usr/bin/env python

## USAGE:
#   normalize rawdata.txt

## ARGUMENTS:
#   rawdata.txt
#   abfile: antibody file to rename antibody names.'ANTIBODY_REFERENCE.csv'
#   -ctrl: not included
## RETURN:
#   3 matrices of normalized data.
#
import argparse
import pandas
import math
import csv
import copy

VERSION="1.0.0"


def supply_args():
    """
    Input arguments
    """
    controls, control_help=define_controls()
    parser = argparse.ArgumentParser(description='Nanostring output RCC readouts (html format) and '
                                                 'converts to raw data tsv')
    parser.add_argument('--rawdata', type=str, help='rawdata.txt')
    parser.add_argument('--abfile', type=argparse.FileType('r'),help='ANTIBODY_REFERENCE.csv')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    #For norm_geomean
    parser.add_argument('-ctrl', action="store_true", default=False, help=control_help)
    args = parser.parse_args()
    return args


VERSION = "2.0.0"


def define_controls():
    """
        Made this function to avoid incorporating global variables
        And this may be more complex in the future.
    Args: no args
    Return: list of controls, and control help string
    """
    controls = ["MCF7", "HCC1954", "BT474", "HeyA8", "MDA468 control", "MDA468control",
                "MDA468+EGF"]

    # make the control help statement for norm geomean ctrl input
    control_help = "Marks all samples to be used as control. Otherwise controls are " + ",".join(
        [str(i) for i in controls])
    return (controls, control_help)

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

# geometric mean
def gm_mean(vector):
    """Get the geometric mean of a vector of numbers
    Args:
        vector
    Returns:
        float
    """
    log_vector = []
    for num in vector:
        log_vector.append(math.log(float(num)))
    num_sum = math.fsum(log_vector)
    g = math.exp(num_sum / len(vector))
    return g


def df2dictarrays(raw_data):
    """Format the raw data input file into usable arrays of Samples, positive and negative
    Args:
        rawdata pandas dataframe
    Returns:
        dict with keys: pos, neg and samples. Values are two dimensional arrays.
    """
    # iterate over the dataframe, adding each line to the appropriate new array.
    samples = []
    pos = []
    neg = []
    pos.append(list(raw_data.columns))
    samples.append(list(raw_data.columns))
    neg.append(list(raw_data.columns))
    for index, row in raw_data.iterrows():
        if "POS" in row.Name:
            pos.append(row.to_list())
        else:
            samples.append(row.to_list())

    # samples are a matrix with rows are antibodies and columns are samples
    output_dict = {"pos": pos, "samples": samples}
    return output_dict


# geomean norm
def geomean_norm(samples, controls, pos, mouseAb, rabbitAb):
    """Normalizes the data
    Args:
        samples: two dimensional samples matrix from get_output
        controls: list of control samples
        pos: two dimensinal pos matrix from get_output
        mouseAb: list of mouse antibodies
        rabbitAb: list of rabbit antibodies
    returns:
        list: list of four two dimensional arrays from after each of the three normalization steps and the log step.
    """
    # samples are matrix of sampls from get_output
    # controls are list of control names (columns) based on the ctrl_flag
    # pos is matrix of positive controls, formatted like the sample matrix
    # mouseAb is a list of mouse antibodies
    # rabbitAb is a list of rabbit antibodies

    norm_factors={}
    sum_of_sums = 0.0
    header = True
    for row in pos:
        if not header:
            count = 0
            num_controls = 0.0
            for col in row:
                if pos[0][count].strip() in controls:
                    sum_of_sums += float(col)
                    num_controls += 1
                count += 1
        header = False
    mean_of_sums = sum_of_sums / num_controls

    for sample_index in range(1, len(pos[0])):
        sample_sum = 0
        for row_index in range(1, len(pos)):
            sample_sum += float(pos[row_index][sample_index])
        cf = mean_of_sums / sample_sum
        norm_factors[pos[0][sample_index]] = [cf]
        for row_index in range(1, len(samples)):
            samples[row_index][sample_index] = float(samples[row_index][sample_index]) * cf

    samples_1 = copy.deepcopy(samples)

    geomeans = []
    for sample_index in range(1, len(samples[0])):
        if samples[0][sample_index].strip() in controls:
            sample_col = []
            for row in samples:
                sample_col.append(row[sample_index])
            sample_col.pop(0) #removes rownames (ab name)
            sample_geomean = gm_mean([x+1 for x in sample_col]) #add one for gm_mean of zeroes
            sample_geomean = sample_geomean-1  # sub one to correct for add one
            geomeans.append(sample_geomean)

    grand_mean = sum(geomeans) / len(geomeans)

    for sample_index in range(1, len(samples[0])):
        sample_col = []
        for row in samples:
            sample_col.append(row[sample_index])
        sample_col.pop(0)
        sample_geomean = gm_mean([x + 1 for x in sample_col])  # add one for gm_mean of zeroes
        sample_geomean = sample_geomean - 1  # sub one to correct for add one
        cf = grand_mean / sample_geomean
        norm_factors[samples[0][sample_index]].append(cf)

        header = True
        for row in samples:
            if header:
                header = False
            else:
                row[sample_index] = row[sample_index] * cf

    samples_2 = copy.deepcopy(samples)

    mouse_igg_index = 0
    rabbit_igg_index = 0
    count = 0
    for row in samples:
        if row[0].strip() == "MmAb-IgG1":
            mouse_igg_index = count
        elif row[0].strip() == "RbAb-IgG":
            rabbit_igg_index = count
        count += 1

    for sample_index in range(1, len(samples[0])):
        mouse_igg = samples[mouse_igg_index][sample_index]
        rabbit_igg = samples[rabbit_igg_index][sample_index]
        #    for row in samples:
        #        if row[0].strip() in mouseAb:
        #            row[sample_index] = row[sample_index] / mouse_igg
        #        elif row[0].strip() in rabbitAb:
        #            row[sample_index] = row[sample_index] / rabbit_igg
        for row in samples:
            if row[0].strip() in mouseAb:
                if row[sample_index] - mouse_igg <= 0:
                    row[sample_index] = 0
                else:
                    row[sample_index] = row[sample_index] - mouse_igg
            elif row[0].strip() in rabbitAb:
                if row[sample_index] - rabbit_igg <= 0:
                    row[sample_index] = 0
                else:
                    row[sample_index] = row[sample_index] - rabbit_igg

    samples_3 = copy.deepcopy(samples)

    for row in samples:
        for index in range(1, len(samples[0])):
            try:
                row[index] = math.log(row[index] + 1, 2)
            except TypeError:
                pass

    output_list = [samples_1, samples_2, samples_3, samples]
    return output_list, norm_factors


def write_norms(norm_dat):
    """Runs the entire norm geomean process from the raw data file to writing the normalized output files.
    Args:
        #norm_dat: list of 4 arrays
    Returns:
        none
    """

    # write the normalized data to a csv file with the name ending in _NORMALIZED.tsv

    normfile_name = "1_ERCC_NORMALIZED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer_1 = csv.writer(csvfile, delimiter="\t")
        norm_writer_1.writerows(norm_dat[0])

    normfile_name = "2_GEOMEAN_NORMALIZED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer_2 = csv.writer(csvfile, delimiter="\t")
        norm_writer_2.writerows(norm_dat[1])

    normfile_name = "3_IGG_SUBTRACTED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer_2 = csv.writer(csvfile, delimiter="\t")
        norm_writer_2.writerows(norm_dat[2])

    normfile_name = "4_LOG_2_NORMALIZED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer = csv.writer(csvfile, delimiter="\t")
        norm_writer.writerows(norm_dat[3])


def main():
    #avoiding global variables for controls and this might get more complicated in the future
    controls, control_help = define_controls()

    args = supply_args()
    # Change long name to something else if necessary
    with args.abfile:
        ab_name, rabbitAb, mouseAb = parse_Ab_ref(args.abfile)

    raw_data = pandas.read_csv(args.rawdata, sep="\t", index_col=False)

    #Prep for dictionary of arrays
    raw_data_short = raw_data.drop(['CodeClass','Accession'], 1)
    split_data = df2dictarrays(raw_data_short)

    # set controls
    if args.ctrl:
        # set controls to all the samples
        controls = split_data["samples"][0][1:]
    else:
        controls

    norm_dat, norm_factors = geomean_norm(samples=split_data["samples"], controls=controls, pos=split_data["pos"], mouseAb=mouseAb, rabbitAb=rabbitAb)
    write_norms(norm_dat)

if __name__ == "__main__":
    main()
