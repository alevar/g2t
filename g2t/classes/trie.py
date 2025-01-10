# build a trie over intervals in a transcriptome

class TrieNode:
    def __init__(self):
        self.children = {}
        self.transcripts = set()  # IDs of transcripts matching this path

class Trie:
    def __init__(self):
        self.root = TrieNode()

    def insert(self, intervals, transcript_id):
        node = self.root
        for interval in intervals:
            if interval not in node.children:
                node.children[interval] = TrieNode()
            node = node.children[interval]
            node.transcripts.add(transcript_id)

    def query(self, chain):
        node = self.root
        for interval in chain:
            if interval not in node.children:
                return set()  # Chain not found
            node = node.children[interval]
        return node.transcripts