#!/bin/python
import os
import pandas as pd
import argparse
import re



def filenames(infile):
    parent = os.path.dirname(infile)
    filestem = os.path.basename(infile).split('.')[0]
    return parent, filestem

def reduce_PO(PO_outfile, reducedPO, numsamples):
    # since iterating over dataframes is very slow; I will reduce the size of the dataframe first...
    # by creating a reduced file containing only lines with num genes = num species

    with open(PO_outfile, 'r') as PO:
        first = PO.readline()
        with open(reducedPO, 'a') as out:
            out.write(first.strip() + '\n')
        for line in PO:
            if line.strip().split('\t')[0] == line.strip().split('\t')[1] == str(numsamples):
                with open(reducedPO, 'a') as out:
                    out.write(line.strip() + '\n')

def flipflopprint(file, pattern) :
    flag=False
    seq = []
    with open(file, 'r') as f:
        for line in f:
            if re.search(rf"^>{pattern}", line):
                flag=True
                continue
            if re.search(rf"^>", line):
                flag=False
            if flag:
                seq.append(line.strip())
    return(''.join(seq))

def PO2fasta(index, row, samplelist, outfiles, peps):
    for header in samplelist:
        outname = os.path.join(outfiles, 'gene_{}.fasta'.format(str(index)))

        with open(outname, 'a') as out:
            pepfile = os.path.join(peps, header.split('.')[0] + '.pep')

            out.write('>{}_{}\n'.format(header.split('.')[0], row[header]))
            #with open(pepfile, 'r') as fas:
            out.write('{}\n'.format(flipflopprint(pepfile, row[header])))

def make_outdir(indir):
    outfiles = o.spath.join(os.path.dirname(indir), POFastas)
    if not os.path.exists(outfiles):
        os.mkdir(outfiles)
    return outfiles

def main():

#### pasre arguments
    parser = argparse.ArgumentParser(description=
            "Filters ProteinOrtho output to only single copy genes in user \
            provided min number of samples and writes corresponding \
            transdecoded protein sequences to fasta files.")
    parser.add_argument('-p', '--PepFiles', type=str, \
            help='directory containing transdecoder *.pep files')
    parser.add_argument('-i', '--POfile', type=str, help='ProteinOrtho \
            output file in *tsv format')
    parser.add_argument('-n', '--numSamples', type=str, \
            help="minimum number of samples required in output")

    args = parser.parse_args()
################
# step 0 define filepart and outnames
    parent, filepart = filenames(args.POfile)
    reduced = os.path.join(parent, filepart + '.reduced.tsv')
# 1st step: reduce pofile,
    reduce_PO(args.POfile, reduced, args.numSamples)
# 2nd: read in reduced file in pandas df

    PO_raw = pd.read_csv(reduced, sep='\t')

    headers = list(PO_raw)[3:]

    # 3rd: make outfile folder for
    outfiles = os.path.join(filenames(args.PepFiles)[0], 'Fastas')

    if not os.path.exists(outfiles):
        os.mkdir(outfiles)

    for index, row in PO_raw.iterrows():
        PO2fasta(index, row, headers, outfiles, args.PepFiles)

    print('done')
### run
if __name__ == '__main__':
	main()
