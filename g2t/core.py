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
        
        # bam reader memoization - these get updated to expedite the graph parsing
        self.cur_bam_seqid = None
        self.cur_node_idx = None

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
                chain = [[x[0],x[1]-1,[tx.get_tid()]] for x in tx.get_chain()] # -1 here to account for the inclusivity rules of the intervals in the transcriptome
                chains.append(chain)

            # partition chains into disjoint sets
            partitioned_chains = partition_chains(chains)
        
            # add partitioned chains to the trie
            self.graphs[ob.get_seqid()].add_from_chains(partitioned_chains)
            
    def convert_read(self, input_record):
        """
        Convert genomic coordinates of a single-end read to transcriptomic coordinates.
        
        Args:
            input_record (pysam.AlignedSegment): Input BAM record to convert.
            
        Returns:
            pysam.AlignedSegment: Converted read in transcriptomic coordinates,
                or None if conversion is not possible.
        
        Raises:
            AssertionError: If transcript lookup or coordinate conversion fails.
        """
        
        if input_record.query_name == "Env_Vpu.13_290_1":
            print("Found read")
        
        # Extract exons from the input record
        exons = extract_exons(input_record)
        if not exons:
            return None
            
        # Update sequence ID tracking if needed
        if input_record.reference_name != self.cur_bam_seqid:
            self.cur_bam_seqid = input_record.reference_name
            self.cur_node_idx = 0
            
        # Find matching transcripts in the graph
        try:
            matching_ids, self.cur_node_idx = self.graphs[input_record.reference_name].find_chain_path(
                exons,
                self.cur_node_idx
            )
        except KeyError:
            return None
        
        # if input_record.query_name == "Env_Vpu.13_290_1":
        #     return None
        # else:
        #     print("Found read")
            
        if not matching_ids:
            return None
            
        # Process each matching transcript
        records = []
        for tid in matching_ids:
            try:
                # Get transcript and compute coordinates
                tx = self.transcriptome.get_by_tid(tid)
                if tx is None:
                    continue
                    
                transcriptomic_start = tx.get_transcriptomic_position(exons[0][0])
                transcriptomic_end = tx.get_transcriptomic_position(exons[-1][1])
                
                if transcriptomic_start is None or transcriptomic_end is None:
                    continue
                    
                # Create new read with transcriptomic coordinates
                record = pysam.AlignedSegment()
                record.query_name = input_record.query_name
                record.query_sequence = input_record.query_sequence
                record.flag = 0
                record.reference_id = self.output_bam_header_map[tid]
                record.reference_start = transcriptomic_start
                record.mapping_quality = 255
                
                # Calculate CIGAR length
                # For a read spanning positions x to y, length should be y - x + 1
                record.cigar = [(0, transcriptomic_end - transcriptomic_start + 1)]
                record.query_qualities = None
                
                records.append(record)
                
            except (KeyError, IndexError) as e:
                continue
                
        return records
            
    def run(self):
        with pysam.AlignmentFile(self.alignment, "rb") as bam:
            try:
                for record in bam.fetch():
                    if record.is_unmapped:
                        continue
                        
                    converted_reads = self.convert_read(record)
                    if converted_reads is None:
                        continue
                    for read in converted_reads:
                        self.output_bam.write(read)

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