import unittest
import sys
import os
from unittest.mock import MagicMock, patch

# Ensure we can import from parent directory
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import ref_utils

class TestIntegration(unittest.TestCase):

    def setUp(self):
        # Mock Entrez if it's not present or even if it is, to avoid network calls
        self.mock_entrez = MagicMock()
        ref_utils.Entrez = self.mock_entrez

    def test_lookup_pubmed_info_success(self):
        # Setup mock behavior
        # esearch returns a handle that supports .close()
        mock_handle = MagicMock()
        self.mock_entrez.esearch.return_value = mock_handle
        self.mock_entrez.esummary.return_value = mock_handle # Re-use or new mock

        
        # We need to mock Entrez.read separately
        # ref_utils calls Entrez.read(search_handle)
        
        # Difficult to mock "handle" reading precisely without staring at code.
        # Let's mock the return value of Entrez.read directly if we can patch it.
        
        with patch('ref_utils.Entrez.read') as mock_read:
            # First call is for esearch
            # Second call is for esummary
            
            mock_read.side_effect = [
                {"IdList": ["12345"]},          # esearch result
                [{                              # esummary result (list of dicts)
                    "Id": "12345",
                    "Title": "Test Article Title",
                    "Source": "Test Journal",
                    "PubDate": "2024 Jan",
                    "AuthorList": ["Tester A"],
                    "DOI": "10.1234/test",
                    "Volume": "1",
                    "Issue": "1",
                    "Pages": "10-20"
                }]
            ]
            
            query = "Test Article Title"
            result = ref_utils.lookup_pubmed_info(query)
            
            self.assertEqual(result['PMID'], "12345")
            self.assertEqual(result['RIS_Title'], "Test Article Title")
            self.assertEqual(result['Year'], 2024)

    def test_lookup_fallback_if_no_match(self):
         with patch('ref_utils.Entrez.read') as mock_read:
            mock_read.side_effect = [
                {"IdList": []} # No IDs found
            ]
            
            query = "ghost paper"
            result = ref_utils.lookup_pubmed_info(query)
            
            self.assertEqual(result['PMID'], "")
            # Should have fallback data if we provided ref_text?
            # Note: lookup_pubmed_info called with just query here, so no detailed fallback expected in result except defaults.

if __name__ == '__main__':
    unittest.main()
