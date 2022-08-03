#!/bin/python
import glob
import os
import pandas as pd
import re
import argparse

# annotation table - Trinity annotation.xls

def prep_table(xls, pref):
    # reduced Trinotate anotation to only non-empty BLAST results
    outfile = os.path.join(os.path.dirname(xls), os.path.basename(xls).split('.')[0] + '.csv')
    sample = pd.read_excel(xls)
    if pref == 'P':
        sample_red = sample[["transcript_id", "prot_coords", "sprot_Top_BLASTP_hit", "gene_ontology_BLASTP"]]
    elif pref == 'X':
        sample_red = sample[["transcript_id", "sprot_Top_BLASTX_hit", "gene_ontology_BLASTX"]]
    sample_red.to_csv(outfile, sep='!', index=True)
    return (outfile)


def annotate_intersection(meme_summary, csv_annot, taxon, pref):
    # takes meme summary of positive results and writes function into the table
    # while writing any results not containing 'taxon' in annotation as putatve
    # contaminants into new outfile

    # outfilenames from infile name and location
    parent = os.path.dirname(meme_summary)
    base = os.path.basename(meme_summary).split('.')[0]
    outmeme = os.path.join(parent, base + '_annotated.csv')
    print(outmeme)
    outcontamination = os.path.join(parent, base + '_putativecontam.csv')
    print(outcontamination)
    checked = []
    cont_list = []

    with open(meme_summary, 'r') as m:
        header = m.readline().strip().split(';')
        new_header = header
        new_header.insert(1, 'function')
        new_header.insert(2, 'GO-terms')
        with open(outmeme, 'w') as out:
            out.write('{}\n'.format(';'.join(new_header)))
            for line in m:
                gene = line.strip()
                if not gene in checked:
                    checked.append(gene)
                    with open(csv_annot, 'r') as anot:
                        header = anot.readline()
                        for line in anot:
                            if gene + '!' in line:
                                if taxon in line:
                                    a = line.strip().split('!')
                                    # blastp = index 4 preferred, if not blastx, if neither to not annotated
                                    if pref == 'P':
                                        if a[3] == '.':
                                            out.write('{};{}\n'.format(gene, 'NA'))
                                            break
                                        else:
                                            symbol = a[3].split('^')[0]
                                            function = a[3].split('=')[1][:-1]
                                            out.write('{};{}_{};{}\n'.format(gene, symbol, function, a[4]))
                                    
                                    elif pref == 'X':
                                        if a[2] == '.' :
                                            out.write('{};{}\n'.format(gene, 'NA'))
                                            break
                                        else:
                                            symbol = a[2].split('^')[0]
                                            function = a[2].split('=')[1]
                                            out.write('{};{}_{};{}\n'.format(gene, symbol, function, a[3]))
                                            #print(symbol, function, numbers)
                                            break
                                else:
                                    cont_list.append(line)

    with open(outcontamination, 'w') as cfile:
        for item in cont_list:
            cfile.write(item)


# test it

#testannot = "C:\\Users\\cpaetzold\\side-projects\\Birthe\\2021_Mar\\Annotation_dnds.xls"
#meme_pos = "C:\\Users\\cpaetzold\\side-projects\\Birthe\\2021_Mar\\positiveintersect.csv"
#meme_neg = "C:\\Users\\cpaetzold\\side-projects\\Birthe\\2021_Mar\\negativeintersect.csv"
def main():
#############
# parse arguments
    parser = argparse.ArgumentParser(description='combines the results of \
            intersect*.py outputs with the functional annotation and the \
            corresponding GO-Terms from the Trinotate output while filtering \
            out any results marked as putative contaminations where the gene \
            did not blast to a user-specified target taxon.')
    parser.add_argument('-a', '--Annotation', type=str, help='Path to \
            Trinotate*.xls annotation result')
    parser.add_argument('-i', '--intersectionFile', type=str,
            help='Path to *.csv file containing intersection of meme or fel\
            results. Default: None', default='None')
    parser.add_argument('-t', '--taxon', type=str, help='Taxon required for \
            annotated gene to be valid. If this taxon is not present in the \
            annotation, the gene will be considered a putative contaminant. \
            Default: Spermatophyta', default='Spermatophyta')
    parser.add_argument('-b', '--preferredBlast', type=str, help='blast result \
            used preferentially for annotation, X or P. Default: P', default='P')

    args = parser.parse_args()

## run it
    print('running prep')
    reduced_annot = prep_table(args.Annotation, args.preferredBlast)
    print('start annot')
    annotate_intersection(args.intersectionFile, reduced_annot, args.taxon, args.preferredBlast)
 

if __name__ == '__main__':
    main()
