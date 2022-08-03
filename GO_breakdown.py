#!/bin/python

from urllib.request import urlretrieve
# Import the OBO parser from GOATools
from goatools import obo_parser
from goatools import base
import os
import re
import argparse
from pathlib import Path

def get_current_obo(directory):
    filename = os.path.join(directory, 'go-basic.obo')
    urlretrieve ("http://geneontology.org/ontology/go-basic.obo", \
                        filename)

    return filename


def go_iterparents(goid, categorydictionary, parsed_obo):
    # pass lowest level goterm to goatools parser (gt = parser)
    gt = parsed_obo[goid]


    while int(gt.level) > 1:

        for term in gt.parents:
            level = str(term.level)
            goterm = term.id
            if level in categorydictionary.keys():
                if goterm in categorydictionary[level].keys():
                    categorydictionary[level][goterm] += 1
                else:
                    categorydictionary[level][goterm] = 1
            else:
                categorydictionary[level] = {}
                categorydictionary[level][goterm] = 1

            gt = parsed_obo[term.id]
    return categorydictionary

def build_dictionary(goterm, level, cat_dict, cat_tracker):

        if level in cat_dict.keys():
            if goterm in cat_dict[level].keys():
                cat_dict[level][goterm] + 1
            else:
                cat_dict[level][goterm] = 1
        else:
            cat_dict[level] = {}
            cat_dict[level][goterm] = 1

        #check level - colect lowest
        if int(level) < cat_tracker[0]:
            cat_tracker = [int(level), goterm]

        return(cat_tracker)

def write_file(outname, cat_dict):
    with open(os.path.join(outname), 'w') as out:
        for level in sorted(cat_dict.keys()):
            out.write('Level {}\n'.format(level))
            for goterm, number in cat_dict[level].items():
                out.write('{};{}\n'.format(goterm, number))


def main():

################################
# argument
    parser = argparse.ArgumentParser(description='summarizes GO-terms per Level \
            for a given annotated intersection, i.e. the results of funct_annotation.py. \
            Produces three tables, one each with all GO-terms per category, \
            biological process, molecular function and cellular compinent.')
    parser.add_argument('-a', '--annotatedFile', type=str, help='Path to \
            *intersect_annotated.csv file')
    parser.add_argument('-o', '--oboFile', type=str, help='Path to go-basic.obo \
            file. If not specified, the current version of the file will be \
            downloaded')
    args = parser.parse_args()





#####################################

    parent = os.path.dirname(args.annotatedFile)
    # download obo database for GO-terms
    if not args.oboFile:
        obofile = get_current_obo(parent)
        go = obo_parser.GODag(obofile)
    else:
        go = obo_parser.GODag(args.oboFile)

    # initialize dictionaries
    mf_dict = {}
    bp_dict = {}
    cc_dict = {}

    with open(args.annotatedFile, 'r') as f:
        for line in f:
            goterms = re.findall('GO:[0-9]*', line)
            if len(goterms) > 0:
                #print(goterms)
                gene = line.strip().split(';')[0]

                # initiate level trackers
                mf = [100, 'xxx']
                bp = [100, 'xxx']
                cc = [100, 'xxx']
                #print(goterms)
                for i in range(len(goterms)):

                    # get highest order term in the list of Trinotate -
                    # then try to get parents of the highest included term
                    # to that end get is level
                    goterm = goterms[i]
                    try:
                        go_id = go[goterms[i]]
                        category = go_id.namespace
                        level = str(go_id.level)
                        if category == 'biological_process':
                            #print('bp')
                            # conditionals to add to specific dictionaries
                            bp = build_dictionary(goterm, level, bp_dict, bp)
                        elif category == 'molecular_function':
                            #print('mf')
                            mf = build_dictionary(goterm, level, mf_dict, mf)
                        elif category == 'cellular_component':
                            #print('cc')
                            cc = build_dictionary(goterm, level, cc_dict, cc)
                    except:
                        # this is for a KeyError, when a GO-term was deprecated
                        # between Uniprot download and the time of running
                        # this script
                        continue

                if not bp[0] == 100:
                    go_iterparents(bp[1], bp_dict, go)

                if not mf[0] == 100:
                    go_iterparents(mf[1], mf_dict, go)

                if not cc[0] == 100:
                    go_iterparents(cc[1], cc_dict, go)

    write_file(os.path.join(parent, Path(args.annotatedFile).stem + '_biolproc.csv'), bp_dict)
    write_file(os.path.join(parent, Path(args.annotatedFile).stem + '_molfunc.csv'), mf_dict)
    write_file(os.path.join(parent, Path(args.annotatedFile).stem + '_celcomp.csv'), cc_dict)

if __name__ == '__main__':
    main()
