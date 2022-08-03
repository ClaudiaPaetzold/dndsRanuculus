#! /bin/python

import os
import glob
from goatools import obo_parser
from goatools import base
import argparse
from pathlib import Path



def main():

################################
# argument
    parser = argparse.ArgumentParser(description='takes the go_breakdown output \
            file as input and adds the decription to each go-term.')
    parser.add_argument('-f', '--inputFile', type=str, help='Path to \
            *biolproc/molfunc/celcomp.csv file')
    parser.add_argument('-o', '--oboFile', type=str, help='Path to go-basic.obo \
            file. If not specified, the current version of the file will be \
            downloaded')
    args = parser.parse_args()





#####################################

    parent = os.path.dirname(args.inputFile)
    outfile = os.path.join(parent, Path(args.inputFile).stem + '_extended.csv')
    # download obo database for GO-terms
    if not args.oboFile:
        obofile = get_current_obo(parent)
        go = obo_parser.GODag(obofile)
    else:
        go = obo_parser.GODag(args.oboFile)


    with open(args.inputFile, 'r') as short:
        with open(outfile, 'w') as long:
            for line in short:
                if 'GO:' in line:
                    goterm, number = line.strip().split(';')
                    go_id = go[goterm]
                    func = go_id.name
                    long.write('{};{};{}\n'.format(goterm, func, number))
                else:
                    long.write(line)


if __name__ == '__main__':
    main()
