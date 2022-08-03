# this script reads in the summary files ceated by summarize_meme.py ...
# ... and creates venn diagrams for the ones with detected positive selection ..
# ... and the ones with where no positive selection was detected

# if you want to use it, replace the filenames and descriptions ...
#  (e.g. hyb1, venn_positve.svg etc) ...
# ... with the approbriate names and filepaths.
# you may also change the colors of the circles, please see here the colors ...
# available : https://matplotlib.org/stable/gallery/color/named_colors.html
# then change e.g. 'tab:red' ot 'darkslateblue' to any color you like


import os
import re
import argparse
from matplotlib import pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn import venn3_circles


def draw_venn3(file1, file2, file3, n1, n2, n3, col1, col2, col3, title, a, outname):
    f1 = []
    with open(file1, 'r') as fn1:
        header = fn1.readline()
        for l1 in fn1:
            try:
                gene = re.search(r'\d+', l1.split(';')[0]).group()
            except:
                gene = l1.split(';')[0]
            f1.append(gene)
        

    f2 = []
    with open(file2, 'r') as fn2:
        header = fn2.readline()
        for l2 in fn2:
            try:
                gene = re.search(r'\d+', l2.split(';')[0]).group()
            except:
                gene = l2.split(';')[0]
            f2.append(gene)
        

    f3 = []
    with open(file3, 'r') as fn3:
        header = fn3.readline()
        for l3 in fn3:
            #print(l3)
            try:
                gene = re.search(r'\d+', l3.split(';')[0]).group()
            except:
                gene = l3.split(';')[0]
            f3.append(gene)
        

    # for the positives

    venn3([set(sorted(f1)), set(sorted(f2)), set(sorted(f3))], (n1, n2, n3), (col1, col2, col3), alpha = a)

    if title:
        plt.title(title, fontsize=20)

    plt.draw()
    plt.savefig('{}.{}'.format(outname, 'svg'), format='svg')
    plt.savefig('{}.{}'.format(outname, 'pdf'), format='pdf')

    #plt.show()


def main():
    #############

        parser = argparse.ArgumentParser(description=
                "Produces Venn diagrams for three samples in *.svg and *.pdf \
                written into the home directory as File1. \
                Please find names of availlable colours here:  \
                https://matplotlib.org/stable/gallery/color/named_colors.html")
        parser.add_argument('--File1', type=str, \
                help="MEME summary file of first sample")
        parser.add_argument('--File2', type=str, \
                help="MEME summary file of second sample")
        parser.add_argument('--File3', type=str, \
                help="MEME summary file of third sample")
        parser.add_argument('--Name1', type=str, \
                help="Name/Title of first sample")
        parser.add_argument('--Name2', type=str, \
                help="Name/Title of second sample")
        parser.add_argument('--Name3', type=str, \
                help="Name/Title of third sample")
        parser.add_argument('--color1', type=str, \
                help="Color for circle corresponding to first sample. \
                Default: darkslateblue", default='darkslateblue')
        parser.add_argument('--color2', type=str, \
                help="Color for circle corresponding to second sample. \
                Default: seagreen", default='seagreen')
        parser.add_argument('--color3', type=str, \
                help="Color for circle corresponding to third sample. \
                Default: hotpink", default='hotpink')
        parser.add_argument('--Title', type=str, \
                help="Title displayed above the diagram. Default: None", \
                default=False)
        parser.add_argument('--transparency', type=int, \
                help="Transparency of circles in diagram. Default: 0.7", \
                default= 0.7)
        parser.add_argument('--outname', type=str, \
                help="Names for outfiles")
        args = parser.parse_args()
    ###################

        basename = os.path.dirname(args.File1)
        outfile = os.path.join(basename, args.outname)
        draw_venn3(args.File1, args.File2, args.File3, args.Name1, args.Name2, \
                    args.Name3, args.color1, args.color2, args.color3, \
                    args.Title, args.transparency, outfile)
if __name__ == '__main__':
    main()
