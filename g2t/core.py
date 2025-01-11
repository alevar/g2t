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
from .classes.graph import Graph

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
        self.transcriptome = Transcriptome()
        self.transcriptome.build_from_file(self.annotation)
        self.transcriptome.coordinate_sort()
        self.transcriptome.build_index()
# 
        self.graphs = dict()
        self.build_interval_graph()

        # initiate output BAM file
        self.output_bam = None
        self.output_bam_header_map = dict()
        self.initiate_output_bam()

    def initiate_output_bam(self):
        # write header using seqids from the transcriptome
        header = {
            'HD': {'VN': '1.6', 'SO': 'coordinate'},
            'SQ': []
        }
        
        # Add each sequence to the header
        for i, tx in enumerate(self.transcriptome):
            header['SQ'].append({
                'SN': tx.get_tid(),
                'LN': tx.elen()
            })
            self.output_bam_header_map[tx.get_tid()] = i

        self.output_bam = pysam.AlignmentFile(self.output, "wb", header=header)
        
    def build_interval_graph(self):
        self.graphs = dict()
        for ob in self.transcriptome.bundle_it(): # iterate over overlap bundles
            if ob.get_seqid() not in self.graphs:
                self.graphs[ob.get_seqid()] = Graph()
            
            # extract chains of overlapping transcripts
            chains = []
            for tx in ob:
                chain = [[x[0],x[1],[tx.get_tid()]] for x in tx.get_chain()]
                chains.append(chain)

            # partition chains into disjoint sets
            partitioned_chains = partition_chains(chains)
        
            # add partitioned chains to the trie
            self.graphs[ob.get_seqid()].add_from_chains(partitioned_chains)

    def run(self):
        with pysam.AlignmentFile(self.alignment, "rb") as bam:
            try:
                for record in bam.fetch():
                    if record.is_unmapped:
                        continue
                        
                    exons = extract_exons(record)
                    if not exons:
                        continue
                    
                    # find matching transcripts in the graph
                    matching_ids = self.graphs[record.reference_name].find_chain_path(exons)
                    if not matching_ids:
                        continue

                    # compute transcriptomic coordinates of the read for each of the matching transcripts
                    for tid in matching_ids:
                        tx = self.transcriptome.get_by_tid(tid)
                        assert tx is not None, f"Transcript {tid} not found in the transcriptome."
                        transcriptomic_start = tx.get_transcriptomic_position(exons[0][0])
                        transcriptomic_end = tx.get_transcriptomic_position(exons[-1][1])
                        assert transcriptomic_start is not None, f"Transcriptomic start position for {tid} not found."
                        assert transcriptomic_end is not None, f"Transcriptomic end position for {tid} not found."

                        # create a new read with the transcriptomic coordinates
                        record = pysam.AlignedSegment()
                        record.query_name = f"{record.query_name}"
                        record.query_sequence = record.query_sequence
                        record.flag = 0
                        record.reference_id = self.output_bam_header_map[tid]
                        record.reference_start = transcriptomic_start
                        record.mapping_quality = 255
                        record.cigar = [(0, transcriptomic_end - transcriptomic_start)]
                        record.next_reference_id = self.output_bam_header_map[tid]
                        record.next_reference_start = transcriptomic_start
                        record.template_length = 0
                        record.query_qualities = None

                        self.output_bam.write(record)

            except Exception as e:
                print(f"Error processing BAM file: {e}")
                return
        
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