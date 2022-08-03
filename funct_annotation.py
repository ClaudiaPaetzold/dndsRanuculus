#!/bin/python
import glob
import os
import pandas as pd
import re
import argparse

# annotation table - Trinity annotation.xls

def prep_table(xls):
    # reduced Trinotate anotation to only non-empty BLAST results
    outfile = os.path.join(os.path.dirname(xls), os.path.basename(xls).split('.')[0] + '.csv')
    sample = pd.read_table(xls, sep='\t')
    sample_red = sample[["transcript_id", "sprot_Top_BLASTX_hit", "prot_coords", "sprot_Top_BLASTP_hit", "gene_ontology_BLASTX", "gene_ontology_BLASTP"]]
    sample_red.to_csv(outfile, sep='!', index=True)
    return (outfile)


def annotate_positives(meme_summary, csv_annot, taxon):
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
        new_header.insert(-1, 'GO-terms')
        with open(outmeme, 'w') as out:
            out.write('{}\n'.format(';'.join(new_header)))
            for line in m:
                memefilename = line.strip().split(';')[0]
                gene = '_'.join(memefilename.split('_')[1:3])
                numbers = line.strip().split(';')[1:]
                if not gene in checked:
                    checked.append(gene)
                    with open(csv_annot, 'r') as anot:
                        header = anot.readline()
                        for line in anot:
                            if gene + '!' in line:
                                if taxon in line:
                                    a = line.strip().split('!')
                                    # blastp = index 4 preferred, if not blastx, if neither to not annotated
                                    if a[4] == '.':
                                        symbol = a[2].split('^')[0]
                                        function = a[2].split('=')[1].split(';')[0]
                                        # go terms blast x index 5
                                        gos = re.findall('GO:[0-9]*', a[5])
                                        out.write('{};{}_{};{};{}\n'.format(gene, symbol, function, ';'.join(numbers), ';'.join(gos)))
                                        break

                                    else:
                                        symbol = a[4].split('^')[0]
                                        function = a[4].split('=')[1].split(';')[0]
                                        # go terms blast x index 6
                                        gos = re.findall('GO:[0-9]*', a[6])
                                        out.write('{};{}_{};{};{}\n'.format(gene, symbol, function, ';'.join(numbers), ';'.join(gos)))
                                        break
                                else:

                                    c = line.strip().split('!')
                                    if c[4] == '.' and c[2] == '.':
                                        out.write('{};{};{}\n'.format(gene, 'NA', ';'.join(numbers)))
                                    else:
                                        cont_list.append(line)

    with open(outcontamination, 'w') as cfile:
        for item in cont_list:
            cfile.write(item)


def annotate_negatives(meme_summary, csv_annot, taxon):
    # same as annotate_positives but other input file (no dnds ratios, thus
    # other output format necessary)

    # outfilenames from infile name and location
    parent = os.path.dirname(meme_summary)
    base = os.path.basename(meme_summary).split('.')[0]
    outmeme = os.path.join(parent, base + '_annotated.csv')
    outcontamination = os.path.join(parent, base + '_putativecontam.csv')


    checked = []
    cont_list = []

    with open(meme_summary, 'r') as m:
        header = m.readline().strip()
        new_header = header
        with open(outmeme, 'w') as out:
            out.write('{};{}\n'.format('gene', 'function'))
            for line in m:
                gene = '_'.join(line.strip().split('_')[1:3])
                #print(gene)
                if not gene in checked:
                    checked.append(gene)
                    with open(csv_annot, 'r') as anot:
                        header = anot.readline()
                        for line in anot:
                            if gene + '$' in line:
                                if taxon in line:
                                    a = line.strip().split('!')
                                    # blastp = index 4 preferred, if not blastx, if neither to not annotated
                                    if a[4] == '.':
                                        if a[2] == '.' :
                                            out.write('{};{}\n'.format(gene, 'NA'))
                                            break
                                        else:
                                            symbol = a[2].split('^')[0]
                                            function = a[2].split('=')[1][:-1]
                                            out.write('{};{}_{}\n'.format(gene, symbol, function))
                                            #print(symbol, function, numbers)
                                            break

                                    else:


                                            symbol = a[4].split('^')[0]

                                            function = a[4].split('=')[1][:-1]
                                            out.write('{};{}_{}\n'.format(gene, symbol, function))
                                            #print(symbol, function, numbers, gene)
                                            break
                                else:

                                    c = line.strip().split('!')
                                    if c[4] == '.' and c[2] == '.':
                                        out.write('{};{}\n'.format(gene, 'NA'))
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
    parser.add_argument('-p', '--positive_intersection', type=str,
            help='Path to *.csv file containing intersection of positive meme \
            results. Default: None', default='None')
    parser.add_argument('-n', '--negative_intersection', type=str, help='Path to *.csv \
            file containing intersection of negative meme results. Default: None',
            default='None')
    parser.add_argument('-t', '--taxon', type=str, help='Taxon required for \
            annotated gene to be valid. If this taxon is not present in the \
            annotation, the gene will be considered a putative contaminant. \
            Default: Spermatophyta', default='Spermatophyta')

    args = parser.parse_args()

## run it
    reduced_annot = prep_table(args.Annotation)
    if not args.positive_intersection == 'None':
        annotate_positives(args.positive_intersection, reduced_annot, args.taxon)
    if not args.negative_intersection == 'None':
        annotate_positives(args.negative_intersection, reduced_annot, args.taxon)

if __name__ == '__main__':
    main()
