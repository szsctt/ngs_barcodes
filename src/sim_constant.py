#!/usr/bin/env python3

# takes as input a barcode yaml file, and simulates reads with
# a random number of combinations of each barcode in each set
# for each set, also incldues a set of reads with different barcode
# (i.e. a barcode that doesn't match the set)
# also specify the total length of reads to simulate

import argparse
import yaml
import random
import itertools

def main():

    parser = argparse.ArgumentParser(description='Simulate reads with five combinations of each barcode in each set')
    parser.add_argument('-b', '--barcode-file', help='barcode yaml file', required=True)
    parser.add_argument('-l', '--read-length', help='total length of reads to simulate', default=150)
    parser.add_argument('-o', '--output-fastq', help='output fastq name', default='simulated_reads.fq')
    parser.add_argument('-c', '--counts', help='output counts for combinations of barcodes', default='simulated_counts.txt')
    args = parser.parse_args()

    # read in barcode yaml file
    with open(args.barcode_file, 'r') as f:
        barcodes = yaml.safe_load(f)

    # check that barcode yaml file is in correct format
    check_barcode_yaml(barcodes)

    # generate reads
    write_reads(barcodes, args.read_length, args.output_fastq, args.counts)
    
def write_reads(barcodes, read_length, output_fastq, output_counts):

    # generate random sequence to fill in between barcodes
    # (length of sequence is read_length - length of barcodes)
    total_length = 0
    for barcode in barcodes:
        bc_set = list(barcode.values())[0]['barcodes']
        total_length += len(list(bc_set.values())[0])

    # to include 'none' barcode, add extra barcode to each set
    for barcode in barcodes:
        bc_set = barcode[list(barcode.keys())[0]]['barcodes']
        bc_set['none'] = generate_random_sequence(len(list(bc_set.values())[0]), 'G')

    # open output file
    with open(output_fastq, 'w') as f, open(output_counts, 'w') as c:
        
        # write header for counts
        header = '\t'.join([list(i.keys())[0] for i in barcodes]) + '\tcount\n'
        c.write(header)

        for bcs, rname, rseq in gen_reads(barcodes, read_length):

            # get a random count for this read
            count = random.randint(1, 100)

            # write to fastq
            for i in range(count):
                f.write(rname + '\n')
                f.write(rseq + '\n')
                f.write('+\n')
                f.write('I'*len(rseq) + '\n')

            # write to counts
            count_line = '\t'.join(bcs) + '\t' + str(count) + '\n'
            c.write(count_line)

def gen_reads(barcodes, read_len):

    parts = []
    pos = 0
    # for each set of barcodes
    for bs in barcodes:
        # get name, start, and barcodes
        bs_name = list(bs.keys())[0]
        bs_start = bs[bs_name]['start']
        
        # add sequence to parts for constant part of barcodes
        parts.append((generate_random_sequence(bs_start - pos),))

        # add barcode names to parts
        parts.append(list(bs[bs_name]['barcodes'].keys()))
       
        # update position in random_sequence
        len_bc = len(list(bs[bs_name]['barcodes'].values())[0])
        pos = bs_start + len_bc 

    # add last part of random sequence
    parts.append((generate_random_sequence(read_len - pos),))

    # generate all combinations of barcodes
    for bc in itertools.product(*parts):
        # get barcode names
        bc_names = [bc[i] for i in range(1, len(bc), 2)]
        # generate read name
        read_name = '@' + '__'.join(bc_names)
        # generate read sequence
        read_parts = []
        for i, part in enumerate(bc):
            if i % 2 == 0:
                read_parts.append(part)
            else:
                read_parts.append(barcodes[i//2][list(barcodes[i//2].keys())[0]]['barcodes'][part])

        yield bc_names, read_name, ''.join(read_parts)


def generate_random_sequence(length, let='A'):
    
    # generate random sequence of length 'length'
    # sequence is composed of 'A' characters only
    return ''.join(let for i in range(length))

def check_barcode_yaml(barcodes):

    # check that barcodes contains a list
    assert isinstance(barcodes, list)

    # check each element of barcodes
    for barcode in barcodes:
        # check that each element of barcodes is a dict
        assert isinstance(barcode, dict)
        # check that length of dict is 1
        assert len(barcode) == 1
        # descend into value of dict
        barcode = barcode[list(barcode.keys())[0]]
        # check that each dict has a 'start' key
        assert 'start' in barcode.keys()
        # check that each dict has a 'barcodes' key
        assert 'barcodes' in barcode.keys()
        # check that each 'barcodes' key contains a dict
        assert isinstance(barcode['barcodes'], dict)
        # check that values of dict are composed of 'ACGT' characters
        for value in barcode['barcodes'].values():
            assert set(value).issubset('ACGT')
        # check that each 'barcodes' key contains at least one barcode
        assert len(barcode['barcodes']) > 0
        # check that all barcode keys are different
        assert len(barcode['barcodes']) == len(set(barcode['barcodes'].keys()))
        # check that all barcode values are different
        assert len(barcode['barcodes']) == len(set(barcode['barcodes'].values()))
        # check that all barcode values are the same length
        assert len(set([len(value) for value in barcode['barcodes'].values()])) == 1


    # check that 'start' for each set are in ascending order
    starts = [list(barcode.values())[0]['start'] for barcode in barcodes]
    assert starts == sorted(starts)


if __name__ == "__main__":
    main()
