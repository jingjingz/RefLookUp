import unittest
import sys
import os
import glob

# Ensure we can import from parent directory (assuming tests/ is one level deep)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import ref_utils

class TestAllFormats(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__), 'test_data')
    
    def test_all_formats_structure(self):
        """
        Iterate over all files in test_data.
        Ensure each parses into a list of dicts with 'title', 'text', 'metadata'.
        Ensure 'text' is a reasonable string.
        """
        files = glob.glob(os.path.join(self.test_dir, '*'))
        if not files:
            self.fail("No test files found in test_data")
            
        for fpath in files:
            with self.subTest(file=os.path.basename(fpath)):
                print(f"Testing {os.path.basename(fpath)}...")
                with open(fpath, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                
                try:
                    entries = ref_utils.parse_any_input(content, fpath)
                except Exception as e:
                    self.fail(f"Exception parsing {os.path.basename(fpath)}: {e}")
                
                self.assertIsInstance(entries, list, f"Output must be a list for {fpath}")
                self.assertTrue(len(entries) > 0, f"Output list empty for {fpath}")
                
                for i, entry in enumerate(entries):
                    # Check standard keys
                    self.assertIn('title', entry, f"Missing 'title' in entry {i} of {fpath}")
                    self.assertIn('text', entry, f"Missing 'text' in entry {i} of {fpath}")
                    self.assertIn('metadata', entry, f"Missing 'metadata' in entry {i} of {fpath}")
                    
                    # Check content types
                    self.assertIsInstance(entry['metadata'], dict, f"'metadata' must be dict in entry {i}")
                    
                    # Specific check for RIS: Text should look like a citation, not raw RIS tag lines
                    if fpath.lower().endswith('.ris'):
                        text = entry['text']
                        # Should not start with "TY  - " or contain raw field tags as the main text
                        # It should look like "Authors. Title..."
                        self.assertFalse(text.strip().startswith("TY  -"), 
                                         f"RIS entry {i} text starts with raw TY tag: {text}")
                        # Check for reasonable length or content
                        self.assertTrue(len(text) > 10, f"RIS entry text too short: {text}")

if __name__ == '__main__':
    unittest.main()
