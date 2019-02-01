import os
import csv

channel_suffix = ['_ext','_p1','_p2','_p3']
        
def getSubDirectoryNames(eosDir):
    # return a list of sub-directory names under eos
    subdirs = os.popen('eos root://cmseos.fnal.gov ls ' + eosDir).read().strip().split('\n')
    return subdirs

def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]


# read dataset list and convert to dictionary
def getDatasetDict(datasetList):
    datadict = {}
    reader = csv.DictReader(open(datasetList, 'rb'))
    for row in reader:
        datadict[row['sample']] = row

    return datadict

# samples to submit crab jobs
def getSamples(infile):
    f=open(infile, 'r')
    tmp=f.read()
    return tmp.replace('\n','').replace(' ','').split(',')


def getSampleFullname(sample, datalist):

    if 'data' in sample:
        sample = sample.replace('fakes_','')
        sample = sample.replace('flips_','')
        sample = sample.replace('obs_','')
        sample = sample.replace('_missingLumis','')
    sample = sample.replace('_jesup','')
    sample = sample.replace('_jesdown','')
    sample = sample.replace('_tesup','')
    sample = sample.replace('_tesdown','')

    if datalist.split('.')[-1]=='csv':
        datalist_dict = getDatasetDict(datalist)
        if sample in datalist_dict:
            dataset = datalist_dict[sample]['dataset']
            return dataset.split('/')[1]
        else:
            raise ValueError("Invalid sample name "+sample)
    else:
        with open(datalist) as f:
            for line in f:
            
                if sample != line.strip():
                    continue

                fullname = f.next().strip()
                return fullname.split('/')[1] # [0] is expected to always be "'"

            print 'WARNING: CANNOT find matched sample name!'
