import os
import re
import copy
import shutil
import argparse
import subprocess

from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

from .utils.common import *

from .classes.txgroup import Transcriptome, Gene, Bundle
from .classes.transcript import Transcript, Object
from .classes.trie import Trie

class G2T:
    def __init__(self, args):
        if not os.path.exists(args.alignment):
            raise FileNotFoundError(f"Alignment file {args.alignment} not found.")
        
        if not os.path.exists(args.annotation):
            raise FileNotFoundError(f"Annotation file {args.annotation} not found.")
            
        if gtf_or_gff(args.annotation) is None:
            raise ValueError(f"{args.annotation} is not a valid GTF/GFF file.")
        
        # INPUT FILES
        # create copies of files in tmp directory for use in the pipeline
        self.alignment = args.alignment
        self.annotation = args.annotation
        self.output = args.output
        self.gtf = gtf_or_gff(args.annotation)
        
        # load transcriptome
        self.transciptome = Transcriptome()
        self.transciptome.build_from_file(self.annotation)
        self.transciptome.coordinate_sort()
        self.transciptome.build_index()
        
        # construct a trie over the transcriptome
        self.trie = self.build_interval_trie()
        
    def build_interval_trie(self):
        trie = Trie()
        for tx in self.transciptome:
            trie.insert(tx.get_chain(), tx.get_tid())

    def run(self):
        print("run")
        
def main():
    parser = argparse.ArgumentParser(description="Tool for HIV-1 genome annotation")

    parser.add_argument('-a', '--annotation', required=True, type=str, help='Path to the reference GTF/GFF annotation file')
    parser.add_argument('-b', '--alignment', required=True, type=str, help='Path to the genomic alignment file to be converted')
    parser.add_argument('-o', '--output', type=str, help='Path to the output BAM file, with mappings in transcriptomic coordinate space')

    args = parser.parse_args()

    g2t = G2T(args)
    g2t.run()

if __name__ == "__main__":
    main()