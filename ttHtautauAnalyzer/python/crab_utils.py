import os 

channel_suffix = ['_ext','_p1','_p2','_p3']

def getSampleFullname(sample, list):

    with open(list) as f:
        for line in f:
            if sample != line.strip():
                continue

            fullname = f.next().strip()
            return fullname.split('/')[1] # [0] is expected to always be "'"


def getSubDirectoryNames(eosDir):
    # return a list of sub-directory names under eos
    subdirs = os.popen('eos root://cmseos.fnal.gov ls ' + eosDir).read().strip().split('\n')
    return subdirs

def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]
