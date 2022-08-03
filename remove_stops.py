import glob
import os
import argparse





def no_stop_codon(file):
    stops = ['TGA', 'TAA', 'TAG']

    parent = os.path.dirname(file)
    outfile = '.'.join(os.path.basename(file).split('.')[:-1]) + '_remove_stop.fasta'
    seqs = {}
    with open(file,'r') as f:
        for line in f:
            if line[0] == '>':
                header = line.strip()
            else :
                sequence = line.strip()
                if len(sequence)%3 != 0:
                    print(file)
                else:
                    seqs[header] = sequence
                    for index in range(0, (len(sequence)), 3):
                        current_codon = sequence[index:index+3]
                        if current_codon in stops:
                            #print(header, index, current_codon)
                            new_seq = sequence[:index] + '---' + sequence[index+3:]
                            if len(new_seq)%3 == 0:
                            #print(new_seq)
                                seqs[header] = new_seq
                            else:
                                print(file, header)
    with open(os.path.join(parent, outfile), 'w') as out:
        for i in seqs:
            out.write('{}\n{}\n'.format(i, seqs[i]))

def main():

###########################
#parse arguments

	parser = argparse.ArgumentParser(description='removes the universal stopcodons from all fasta files in the specified directory')
	parser.add_argument('Path', type=str, help='directory path to fasta files - these should be named *.fasta')

	args = parser.parse_args()


###########################
	for file in glob.glob(os.path.join(args.Path, '*.fasta')):
		no_stop_codon(file)

if __name__ == '__main__':
	main()
