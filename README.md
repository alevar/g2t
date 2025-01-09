# g2t
Genome To Transcriptome re-mapping utility

Performs conversion of genomic coordinates in alignment files into transcriptomic. Augments mappings with all transcriptomic multimappers. For each genomic chain - searches an Interval Trie indexed with an Array-Backed Interval Tree for a list of transcripts containing full chain.

Currently implemented in Python, with plans to transition to Rust once stable
