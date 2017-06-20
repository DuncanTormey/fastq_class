#!/usr/bin/env python
#Author: Duncan Tormey
#Email: dut@stowers.org or duncantormey@gmail.com

from __future__ import print_function
import os
import re
import glob
import pandas as pd
from pprint import pprint

def read_list_file(path):
    l = list(set(open(path.strip('\n'), 'r').read().strip('\n').splitlines()))
    return l

class Fastq(object):
    def __init__(self, fastq_path):
        self.path = fastq_path
        self.file_name = os.path.basename(self.path)
        if re.match('s_[0-9]+_[0-9]+_', self.file_name[:6]):
            self.lane = self.file_name[2]
            self.pair = self.file_name[4]

    def __repr__(self):
        return self.path


class SimrFastq(Fastq):
    def __init__(self, fastq_path, sample_report_path):
        super(SimrFastq, self).__init__(fastq_path)
        self.flowcell = sample_report_path.split('/')[-2]
        self.sample_report_path = sample_report_path
        self.sr_df = pd.read_csv(self.sample_report_path)
        self.sr_df.columns = [c.replace(' ', '_') for c in self.sr_df]
        for column in self.sr_df.columns:
            try:
                setattr(self, column,
                        self.sr_df.loc[self.sr_df['output'] == self.file_name][column].values[0])
            except IndexError:
                print('%s %s missing' % (self.path, str(column)))


class SimrData(object):
    '''
    This class takes as input a list of sample names and a list of molng directory paths 
    '''
    def __init__(self, sample_name_list, molng_path_list):
        self.molng_path_list = molng_path_list
        self.sample_name_list =sample_name_list
        self.flowcells = []
        self.fastqs = []
        for path in molng_path_list:
            flowcells = os.listdir(path)
            self.flowcells.extend(flowcells)
            for flowcell in flowcells:
                sample_report_path = '/'.join([path, flowcell, 'Sample_Report.csv'])
                fastq_paths = glob.glob(path + '/' + flowcell + '/*.fastq.gz')
                fastq_paths = [fq for fq in fastq_paths if 'Undetermined' not in fq]
                fastqs = [SimrFastq(fq_path, sample_report_path) for fq_path in fastq_paths ]
                fastqs = [fastq for fastq in fastqs
                          if fastq.sample_name in self.sample_name_list]
                self.fastqs.extend(fastqs)
        self.sample_fastqs = {sample: [fastq for fastq in self.fastqs if fastq.sample_name == sample]
                              for sample in self.sample_name_list}

    def select_flowcells(self,flowcells):
        print('selecting data from flowcells:')
        pprint(flowcells)
        for sample in self.sample_fastqs:
            self.sample_fastqs[sample] = [fq for fq in self.sample_fastqs[sample] if fq.flowcell in flowcells]


    def ret_seq_direct_samples(self):
        sds_sample_fastqs = {}
        for sample in self.sample_fastqs:
            try:
                left_reads = [fq.path for fq in self.sample_fastqs[sample] if fq.read == 1]
                right_reads = [fq.path for fq in self.sample_fastqs[sample] if fq.read  == 2]
            except AttributeError:
                print('Unable to determine sequencing direction from sample report data')

            sds_sample_fastqs[sample] = {'left_reads': sorted(left_reads), 'right_reads': sorted(right_reads)}

        return sds_sample_fastqs

    def ret_seq_direct_merged(self):
        sdm_fastqs = {'left_reads': [], 'right_reads': []}
        for sample in self.sample_fastqs:
            sdm_fastqs['left_reads'].extend([fq.path for fq in self.sample_fastqs[sample] if fq.read == 1])
            sdm_fastqs['right_reads'].extend([fq.path for fq in self.sample_fastqs[sample] if fq.read == 2])

        sdm_fastqs['left_reads'] = sorted(sdm_fastqs['left_reads'])
        sdm_fastqs['right_reads'] =  sorted(sdm_fastqs['right_reads'])

        return sdm_fastqs

if __name__ == "__main__":
    molngs = ['/n/analysis/Baumann/lpn/MOLNG-1608/']
    samples = ['A_R1', 'A_R2', 'A_R3']
    test = SimrData(samples, molngs)
    pprint(test.sample_fastqs)

