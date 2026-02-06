import unittest
import sys
import os

# Ensure we can import from parent directory
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import ref_utils

class TestFormatting(unittest.TestCase):
    
    def test_extract_year(self):
        self.assertEqual(ref_utils.extract_year("Journal of Stuff. 2024;12(1):100"), 2024)
        self.assertEqual(ref_utils.extract_year("No year here"), 0)
        self.assertEqual(ref_utils.extract_year("Published in 1999."), 1999)

    def test_extract_title(self):
        ref = "Smith J. Some Great Paper Title. Journal of Science. 2024."
        self.assertEqual(ref_utils.extract_title(ref), "Some Great Paper Title")
        
        ref_sparse = "Just A Title"
        self.assertEqual(ref_utils.extract_title(ref_sparse), "Just A Title")

    def test_format_authors_nlm(self):
        authors = ["Smith James", "Doe, John"]
        # The function format_authors assumes cleaned list or handles it?
        # Let's check implementation: if style="NLM", it just joins with comma.
        # But wait, logic in format_bibtex_authors does the heavy lifting for formatting names usually.
        # format_authors just joins.
        self.assertEqual(ref_utils.format_authors(authors, "NLM"), "Smith James, Doe, John")

    def test_format_authors_apa(self):
        # "Smith James" -> "Smith, J." logic is inside format_authors for APA
        authors = ["Smith James", "Doe John"]
        expected = "Smith, J. & Doe, J."
        self.assertEqual(ref_utils.format_authors(authors, "APA"), expected)

if __name__ == '__main__':
    unittest.main()
