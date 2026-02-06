import unittest
import sys
import os
import glob

# Ensure we can import from parent directory
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import ref_utils

class TestParsing(unittest.TestCase):

    def test_parse_references_text_simple(self):
        text = "Ref 1\nRef 2\n"
        entries = ref_utils.parse_references_text(text)
        self.assertEqual(len(entries), 2)
        self.assertEqual(entries[0]['text'], "Ref 1")
    
    def test_parse_csv_content(self):
        csv_text = "Title,Authors,Year\nMy Title,Me,2024"
        entries = ref_utils.parse_csv_content(csv_text)
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0]['title'], "My Title")
        self.assertEqual(entries[0]['metadata']['Year'], "2024")

    def test_bibtex_parsing(self):
        bib_text = """
@article{key1,
  title={My Article},
  author={Smith, John},
  journal={Journal of Testing},
  year={2023}
}
"""
        entries = ref_utils.parse_bibtex_content(bib_text)
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0]['title'], "My Article")
        # Check author formatting logic
        self.assertIn("Smith J", entries[0]['metadata']['Authors'])

    def test_real_files_if_exist(self):
        # This test iterates over files in test_data and ensures they parse without error
        test_dir = os.path.join(os.path.dirname(__file__), 'test_data')
        files = glob.glob(os.path.join(test_dir, '*'))
        
        if not files:
            print("\nAggregated test: No files found in test_data, skipping.")
            return

        for fpath in files:
            with self.subTest(file=os.path.basename(fpath)):
                print(f"Testing parsing of {os.path.basename(fpath)}...")
                with open(fpath, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                
                try:
                    entries = ref_utils.parse_any_input(content, fpath)
                    self.assertIsInstance(entries, list)
                    # We expect at least one entry if file is not empty/malformed
                    if content.strip():
                        self.assertGreater(len(entries), 0, f"No entries parsed from {fpath}")
                        
                        # Basic check on first entry
                        first = entries[0]
                        self.assertIn('title', first)
                        self.assertIn('text', first)
                        
                except Exception as e:
                    self.fail(f"Parsing failed for {fpath}: {e}")

    def test_ris_parsing(self):
        ris_text = """
TY  - JOUR
TI  - RIS Title
AU  - Author, One
PY  - 2024
ER  - 
"""
        entries = ref_utils.parse_ris_content(ris_text)
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0]['title'], "RIS Title")
        self.assertIn("Author, One", entries[0]['text'])
        self.assertEqual(entries[0]['metadata']['Year'], "2024")

if __name__ == '__main__':
    unittest.main()
