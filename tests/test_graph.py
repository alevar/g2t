# testing the graph API
# tests/test_your_class.py

import unittest
from g2t.classes.graph import Graph

class Test_Graph(unittest.TestCase):
    def test_greet(self):
        # construct a trie over the transcriptome
        self.graph = Graph()
        chains = [(0, 30, {1,2}), (31, 50, {1}), (60, 80, {1,2}),(90,100, {2})]
        self.graph.add_from_chains(chains)

        # Query path
        path = [[10, 30],[60, 70]]
        matching_ids = self.graph.find_chain_path(path)
        print(matching_ids == {2},matching_ids)

        path = [[65, 75]]
        matching_ids = self.graph.find_chain_path(path)
        print(matching_ids == {1,2},matching_ids)

        path = [[60, 80],[90, 100]]
        matching_ids = self.graph.find_chain_path(path)
        print(matching_ids == {2},matching_ids)
        
        path = [[10, 50],[60, 70]]
        matching_ids = self.graph.find_chain_path(path)
        print(matching_ids == {1},matching_ids)

    def setUp(self):
        self.graph = Graph()

    def test_empty_graph(self):
        """Test behavior with empty graph"""
        target_chain = [(1, 5, {1})]
        self.assertEqual(self.graph.find_chain_path(target_chain), set())

    def test_empty_target(self):
        """Test behavior with empty target chain"""
        self.graph.add_from_chains([(1, 5, {1, 2})])
        self.assertEqual(self.graph.find_chain_path([]), set())

    def test_single_interval_match(self):
        """Test finding a single interval that exactly matches"""
        chains = [
            (1, 5, {1, 2}),
            (6, 10, {2, 3})
        ]
        self.graph.add_from_chains(chains)
        target = [(1, 5)]
        self.assertEqual(self.graph.find_chain_path(target)[0], {1, 2})

    def test_split_interval(self):
        """Test finding a path that requires combining multiple intervals"""
        chains = [
            (1, 3, {1, 2}),
            (4, 6, {1, 2}),
            (7, 10, {2, 3})
        ]
        self.graph.add_from_chains(chains)
        target = [(1, 6)]
        self.assertEqual(self.graph.find_chain_path(target)[0], {1, 2})

    def test_overlapping_intervals(self):
        """Test handling overlapping intervals in the graph"""
        chains = [
            (1, 4, {1, 2}),
            (3, 5, {3}),
            (5, 8, {2, 3})
        ]
        self.graph.add_from_chains(chains)
        target = [(1, 6)]
        self.assertEqual(self.graph.find_chain_path(target)[0], {2})

    def test_id_intersection(self):
        """Test proper ID set intersection along the path"""
        chains = [
            (1, 3, {1, 2, 3}),
            (4, 6, {2, 3, 4}),
            (7, 9, {3, 4, 5})
        ]
        self.graph.add_from_chains(chains)
        target = [(1, 6)]
        self.assertEqual(self.graph.find_chain_path(target)[0], {2, 3})

    def test_no_valid_path(self):
        """Test when no valid path exists"""
        chains = [
            (1, 3, {1, 2}),
            (5, 7, {1, 2}),  # Gap between 3 and 5
            (8, 10, {2, 3})
        ]
        self.graph.add_from_chains(chains)
        target = [(1, 7)]
        self.assertEqual(self.graph.find_chain_path(target)[0], set())

if __name__ == "__main__":
    unittest.main()
