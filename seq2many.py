#!/usr/bin/python3
"""
MIT License

Copyright (c) 2022 Shinji Iida

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import argparse
from argparse import RawTextHelpFormatter
import sys

def ref2mut(ref_seq, position, aa):
    """
    Args: 
        ref_seq: a reference sequence (str)
        position: position (0-index) in the reference sequence that mutate (int)
        aa: one-letter amino acid into which the position mutates (str)

    Return: 
        seq: a mutated sequence (str)
    """
    # NOTE: the reference sequence is duplicated if the specified amino acid is the one in the reference.
    seq = ref_seq[:position] + aa.upper() + ref_seq[position+1:]
    return seq

def read_seq(input_seq):
    ref_seq = ""
    with open(input_seq, "r") as fin:
        for line in fin:
            if line.lstrip().startswith(">"):
                continue
            ref_seq = ref_seq + line.strip()
    return ref_seq 

def write_seq(seqs):

    # Single sequence is outputed.
    if isinstance(seqs, str):
        with open("out.seq", "w") as fout:
            fout.write(f"> mutant \n")
            fout.write(seqs)

    # Many mutant sequences are outputed.
    elif isinstance(seqs, list): 
        for i, seq in enumerate(seqs):
            with open(f"mutant_{i}.seq", "w") as fout:
                fout.write(f"> mutant {i}\n")
                fout.write(seq+"\n")

def output_position_index(seq, prefix):
    with open(f"index_{prefix}.out", "w") as fout:
        fout.write("#index, aa\n")
        for i, seq in enumerate(seq):
            fout.write(f"{i} {seq}\n")

def multiple_mutation(seq, aa_position):
    """
    Args:
        seq (str): a sequence to mutate
        aa_position (dict): pairs between a one-letter amino acid and a 0-indexed position,
            e.g. {'A':2, 'L':10} 
            In this case, the amino acids in seq[2] and seq[10] mutate into A and L.

    Return:
        multseq (str): the mutated sequence 
    """
    if not aa_position: 
        sys.exit(f"ERROR: No items in the dictionary")

    multseq = seq
    for aa, position in aa_position.items():
        multseq = ref2mut(multseq, position, aa)
    return multseq

def deep_mutational_scanning(seq, begin_position, end_position):
    AAs = list("G P A V L I M C F Y W H K R Q N E D S T".split())
    dms_seq = []
    icou = 0
    for i in range(begin_position, end_position):
        for aa in AAs:
            seq = ref2mut(seq, i, aa)
            output_position_index(seq, f"mutant_{icou}")
            dms_seq.append(seq)
            icou += 1
    return dms_seq

def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2), f"The length of each sequqnce is different: {len(seq1)} and {len(seq2)}"
    hdist = 0
    for a1, a2 in zip(seq1, seq2):
        if a1 != a2:
            hdist += 1
    return hdist

def read_mutation_list(mutation_list):
    """
    Args:
        mutation_list: A control file that operates which amino acid in a position mutates into another. 
        
    Return:
        aa_position (dict):  e.g. {'A':2, 'L':10}
    """
    with open(mutation_list, "r") as fin:
        aa_position = {}
        for line in fin:
            aa, position = line.strip().split(":")
            aa_position[aa] = int(position)
    return aa_position

def _len_eq_1(char):
    if len(char) == 1:
        return char
    raise argparse.ArgumentTypeError(f"The input value must be a single character: {char}")

def main():
    p = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    p.add_argument("-i", "--ref", required=True)

    p.add_argument("-m", "--mode", default="single", help="single/multi/deep")
    p.add_argument("-p", "--position", type=int, help="0-index position in the reference sequence")

    grp2 = p.add_mutually_exclusive_group()
    grp2.add_argument("-ia", "--aa", type=_len_eq_1, help="single-letter amino acid")
    grp2.add_argument("-ml", "--mlist", help="Specify a mutation list file: e.g. The file should be formatted like \n A:1\n L:2\n K:5\nThe first line indicates that the first index's amino acid mutates into Alanine, etc")
    grp2.add_argument("-rg", "--region", type=int, action='append', nargs=2, help="Regions (begin and end of 0-index residue No.)")

    args = p.parse_args()
    ref = args.ref
    mode = args.mode
    position = args.position
    aa = args.aa
    mutation_list = args.mlist
    regions = args.region

    ref_seq = read_seq(ref)
    output_position_index(ref_seq, "ref")

    if mode == "single":
        seq = ref2mut(ref_seq, position, aa)
        write_seq(seq)
        output_position_index(seq, mode)

    elif mode == "multi":
        aa_position = read_mutation_list(mutation_list)
        seq = multiple_mutation(ref_seq, aa_position)
        write_seq(seq)
        output_position_index(seq, mode)

    elif mode == "deep":
        for region in regions:
            begin, end = region
            seq = deep_mutational_scanning(ref_seq, begin, end)
            write_seq(seq)

    else:
        print("Invalid mode was specified")

if __name__ == "__main__":
    main()
