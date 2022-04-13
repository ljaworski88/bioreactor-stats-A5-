#!/usr/bin/env python3

import csv
import re
import numpy as np
import pandas as pd
import argparse

def rna_quantization(filename):
    '''
    This small function takes the columnwise data produced by the ___
    microplate reader using the SpectraDrop RNA DNA abs protocol. It then uses
    the sample names to produce a dataframe and xlsx file containing all raw
    concentration values, A260: A280 and A230 ratios, and the averaged
    concentration along with purity reports.
    '''
    with open(filename, 'r', encoding='utf-16-le') as RNA_data:
        RNA_reader = csv.reader(RNA_data, dialect='excel-tab')
        # skipping the plate layout portion of the datafile, it doesn't have
        # what we want
        for row in RNA_reader:
            if row:
                if row[0] == '~End':
                    break
        # making the ording of the exprimental conditions make sense in the output
        day_ordering = pd.api.types.CategoricalDtype(categories=['d0', 'd1', 'd7', 'd14'],
                                                     ordered=True)
        oxygen_ordering = pd.api.types.CategoricalDtype(categories=['x', '5', '21'],
                                                     ordered=True)
        load_ordering = pd.api.types.CategoricalDtype(categories=['x', 'l'],
                                                     ordered=True)
        rna_data = pd.DataFrame(columns=['pig',
                                         'day',
                                         'oxygen',
                                         'loading',
                                         'tissue',
                                         'concentration',
                                         'A260:A280',
                                         'A260:A230',
                                         'avg_concentration'])

        # we have to set the categories before setting the index, it's just easier this way
        rna_data['day'] = rna_data.day.astype(day_ordering)
        rna_data['oxygen'] = rna_data.oxygen.astype(oxygen_ordering)
        rna_data['loading'] = rna_data.loading.astype(load_ordering)

        rna_data.set_index(['pig','day','oxygen','loading','tissue'], inplace=True)

        # The regex for separating the samples by id and groups
        sample_id_regex = re.compile(r'(p[0-9][0-9])-(d[0-9]+)-(x|5|21)-(x|l)-(af|np)')
        sample_id = ()
        # we can now go through the data portion of the file
        for row in RNA_reader:
            # all the rows that we want have more than 7 columans, so we can skip the others
            if len(row) > 7:
                # the sample name is the first column and conforms to the regex
                if sample_id_regex.match(row[0]):
                    sample_id = sample_id_regex.match(row[0])
                    # using the sample name we can split it into columns for use
                    # as an identifier, because every combination of condtions
                    # and tissues is unique. The actual data are contatained in
                    # 'row' and are put into numpy arrays for later processing.
                    data_point = pd.DataFrame({'pig' : [sample_id.group(1)],
                                               'day' : [sample_id.group(2)],
                                               'oxygen' : [sample_id.group(3)],
                                               'loading' : [sample_id.group(4)],
                                               'tissue' : [sample_id.group(5)],
                                               'concentration' : [np.array(row[7], dtype=np.float32)],
                                               'A260:A280' : [np.array(row[5], dtype=np.float32)],
                                               'A260:A230' : [np.array(row[6], dtype=np.float32)]})
                    data_point.set_index(['pig','day','oxygen','loading','tissue'], inplace=True)
                    rna_data = rna_data.append(data_point, sort=True)
                else:
                    # Repeat measurements do not list the sample name again
                    if sample_id:
                        rna_data.loc[sample_id.group(1, 2, 3, 4, 5), 'concentration'] = np.append(rna_data.loc[sample_id.group(1, 2, 3, 4, 5), 'concentration'], np.array(row[7], dtype=np.float32))
                        rna_data.loc[sample_id.group(1, 2, 3, 4, 5), 'A260:A280'] = np.append(rna_data.loc[sample_id.group(1, 2, 3, 4, 5), 'A260:A280'], np.array(row[5], dtype=np.float32))
                        rna_data.loc[sample_id.group(1, 2, 3, 4, 5), 'A260:A230'] = np.append(rna_data.loc[sample_id.group(1, 2, 3, 4, 5), 'A260:A230'], np.array(row[6], dtype=np.float32))
    # get the average concentration of the measurements and check the purity of
    # the measurements
    rna_data['avg_concentration'] = rna_data.concentration.apply(np.mean)
    rna_data['A280-contaminated'] = rna_data['A260:A280'].apply(lambda x: np.any(x < 1.2))
    rna_data['A230-contaminated'] = rna_data['A260:A230'].apply(lambda x: np.any(x < 1.2))
    # when making cDNA we will be using 200ng/uL @ 10uL to load so we will
    # calculate how much of each sample I need in each reaction and how much water
    rna_data['uL_rna'] = 10 * 200  / rna_data.avg_concentration
    rna_data['uL_water'] = 10 - rna_data.uL_rna
    rna_data = rna_data.astype({'A230-contaminated': 'category',
                                'A280-contaminated': 'category',
                                'avg_concentration': 'float64',
                                'uL_rna': 'float64',
                                'uL_water': 'float64'})
    # round to 2 decimal places because the pipettor can't do better than that
    rna_data = rna_data.round({"uL_rna": 2, "uL_water": 2})
    return rna_data

def main():
    # add a cli for the program
    parser = argparse.ArgumentParser(description='quantize RNA reading results')
    parser.add_argument('input', nargs='*')
    parser.add_argument('-o', '--output', nargs='?', default='result.xlsx')
    parser.add_argument('-l', '--loading')
    cli_args = parser.parse_args()

    if cli_args.output is None:
        cli_args.output = 'result.xlsx'
    rna_analysis = pd.DataFrame()
    print('Input File(s):')
    for datafile in cli_args.input:
        if re.search(r'\.txt', datafile):
            print(datafile)
            rna_analysis = rna_analysis.append(rna_quantization(datafile))
        rna_analysis.sort_index(inplace=True)

    print('\nOuput File(s):')
    if re.search(r'\.xlsx', cli_args.output):
        rna_analysis.to_excel(cli_args.output)
        print(cli_args.output)
    else:
        print('-l opition must be a .xlsx file')

    if cli_args.loading is not None:
        if re.search(r'\.xlsx', cli_args.loading):
            rna_analysis.loc[:, ['uL_rna','uL_water']].to_excel(cli_args.loading)
            print(cli_args.loading)
        else:
            print('-l opition must be a .xlsx file')

if __name__ == '__main__':
    main()
