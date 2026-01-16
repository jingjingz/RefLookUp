import re
import csv
import time
import ssl
from datetime import datetime

try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context

try:
    from Bio import Entrez
except ImportError:
    Entrez = None


def parse_references(text):
    """
    Split the input text into references.
    """
    entries = [line.strip() for line in text.splitlines() if line.strip()]
    return entries

def extract_year(reference):
    """Find the first 4‑digit year in the reference; return as int or 0."""
    match = re.search(r'\b(19|20)\d{2}\b', reference)
    return int(match.group()) if match else 0

def extract_title(reference):
    """
    Extract a plausible title from the reference.
    Heuristic: Standard key citations are 'Authors. Title. Journal...'
    We split by '. ' and assume the second part is the title.
    """
    parts = reference.split('. ')
    if len(parts) >= 3:
        title = parts[1].strip()
        if title.endswith('.'):
             title = title[:-1]
        return title
    
    if len(parts) == 2:
        return parts[0].strip()

    return reference.strip()

def format_authors(author_list, style="NLM"):
    """
    Format author list based on style.
    NLM: Smith J, Doe A
    APA: Smith, J., & Doe, A. (if < 20 authors usually)
    """
    if not author_list:
        return ""
    
    if style == "NLM":
        return ", ".join(author_list)
    
    if style == "APA":
        # Entrez usually gives "Smith J". APA needs "Smith, J."
        formatted = []
        for auth in author_list:
            parts = auth.split(" ")
            if len(parts) >= 2:
                last = parts[0]
                # Initials are the rest
                initials = "".join([f"{p[0]}." for p in parts[1:]])
                formatted.append(f"{last}, {initials}")
            else:
                formatted.append(auth)
        
        if len(formatted) == 1:
            return formatted[0]
        elif len(formatted) == 2:
            return f"{formatted[0]} & {formatted[1]}"
        else:
            # APA 7th: list up to 20 authors. For simplicity, we list all or truncate with ... if huge
            return ", ".join(formatted[:-1]) + ", & " + formatted[-1]

    return ", ".join(author_list)

def format_nlm_citation(item):
    """
    Format metadata into NLM style:
    Authors. Title. Journal. Date;Vol(Issue):Pages.
    """
    # Authors
    authors = item.get("AuthorList", [])
    if not authors and "LastAuthor" in item:
        authors = [item["LastAuthor"]]
    
    auth_str = format_authors(authors, "NLM")
    if auth_str:
        auth_str += "."
    
    # Title
    title = item.get("Title", "")
    if title and not title.endswith("."):
        title += "."
    
    # Journal
    journal = item.get("Source", "") or item.get("FullJournalName", "")
    if journal and not journal.endswith("."):
        journal += "."
        
    # Date/Vol/Issue/Pages
    date = item.get("PubDate", "")
    vol = item.get("Volume", "")
    issue = item.get("Issue", "")
    pages = item.get("Pages", "")
    
    # NLM Format: Date;Vol(Issue):Pages.
    # e.g. 2025 Nov 6;80(12):glaf206.
    
    pub_part = ""
    if date:
        pub_part += f"{date}"
    if vol:
        pub_part += f";{vol}"
    if issue:
        pub_part += f"({issue})"
    if pages:
        pub_part += f":{pages}"
    
    if pub_part and not pub_part.endswith("."):
        pub_part += "."
        
    components = [auth_str, title, journal, pub_part]
    return " ".join([c for c in components if c.strip()])

def format_apa_citation(item):
    """
    Format metadata into APA style (approximate 7th ed):
    Authors. (Year). Title. Journal, Vol(Issue), Pages. DOI
    """
    # Authors
    authors = item.get("AuthorList", [])
    if not authors and "LastAuthor" in item:
        authors = [item["LastAuthor"]]
        
    auth_str = format_authors(authors, "APA")
    
    # Year
    pub_date = item.get("PubDate", "")
    # Try extract year
    year_match = re.search(r'\b(19|20)\d{2}\b', pub_date)
    year = year_match.group() if year_match else "n.d."
    
    # Title
    title = item.get("Title", "")
    if title and not title.endswith("."):
        title += "."
        
    # Journal
    journal = item.get("Source", "") or item.get("FullJournalName", "")
    
    # Volume/Issue/Pages
    vol = item.get("Volume", "")
    issue = item.get("Issue", "")
    pages = item.get("Pages", "")
    
    # Construct: Auth. (Year). Title. Journal, Vol(Issue), Pages.
    cit = ""
    if auth_str:
        cit += f"{auth_str} "
    
    cit += f"({year}). {title} "
    
    if journal:
        cit += f"{journal}"
        
    if vol:
        cit += f", {vol}"
        if issue:
            cit += f"({issue})"
            
    if pages:
        cit += f", {pages}"
        
    if not cit.endswith("."):
        cit += "."
        
    return cit

def get_category(item):
    """
    Determine category: Article, Abstract, or Other.
    """
    ptypes = item.get("PubTypeList", [])
    ptypes_lower = [p.lower() for p in ptypes]
    if any("abstract" in p for p in ptypes_lower):
        return "Abstracts/Posters"
    return "Articles" 

def lookup_pubmed_info(query, full_ref_text=None):
    """
    Query PubMed. Returns a dictionary with detailed metadata and IDs.
    Validates result against query title and (optionally) full_ref_text for authors.
    """
    default_res = {
        "OriginalQuery": query,
        "RefText_NLM": "",
        "RefText_APA": "",
        "PMID": "", "PMCID": "", "DOI": "", "Link": "",
        "Year": 0,
        "Category": "Uncategorized"
    }

    if not Entrez:
        return default_res

    Entrez.email = "A.N.Other@example.com"
    
    try:
        # 1. Search - fetch candidates
        time.sleep(0.35)
        # Fetch a few to check for best match
        search_handle = Entrez.esearch(db="pubmed", term=query, retmax=3)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_results["IdList"]
        if not id_list:
            return default_res
            
        # 2. Fetch Summaries for candidates
        summary_handle = Entrez.esummary(db="pubmed", id=",".join(id_list))
        summary_results = Entrez.read(summary_handle)
        summary_handle.close()
        
        if not summary_results:
            return default_res

        from difflib import SequenceMatcher

        best_match = None
        best_score = 0.0

        clean_query = re.sub(r'[^a-zA-Z0-9]', '', query).lower()

        for data in summary_results:
            found_title = data.get("Title", "")
            clean_found = re.sub(r'[^a-zA-Z0-9]', '', found_title).lower()
            
            # 1. Strict Title Check
            # Ratio of 0.9 allows for minor typos/spacing but is essentially "Exact"
            similarity = SequenceMatcher(None, clean_found, clean_query).ratio()
            
            # Also check substring match (sometimes query is cut off or data has extra info)
            if clean_query in clean_found or clean_found in clean_query:
                # Boost if substring match exists
                if similarity < 0.9: 
                    similarity = 0.85 # Artificial boost if it's a substring match
            
            if similarity < 0.85: # Strict threshold
                continue
                
            # 2. Author Check (if reference text available)
            # This handles "same title, different paper"
            if full_ref_text:
                authors = data.get("AuthorList", [])
                if not authors and "LastAuthor" in data:
                    authors = [data["LastAuthor"]]
                
                # Check if at least one author's surname is in the full text
                # We do a loose check: "Zou" in "Zou J, ..."
                author_match = False
                if not authors:
                     # No authors in metadata? Accept if title is very good
                     author_match = True 
                else:
                    for auth in authors:
                        # Extract surname (simple split)
                        surname = auth.split(" ")[0]
                        # Clean surname
                        surname = re.sub(r'[^a-zA-Z]', '', surname)
                        if len(surname) > 1 and surname in full_ref_text:
                            author_match = True
                            break
                
                if not author_match:
                    # Title matched but authors didn't. Likely wrong paper.
                    continue

            # If we passed checks, is this the best one?
            if similarity > best_score:
                best_score = similarity
                best_match = data
        
        if not best_match:
            print(f"  No valid match for '{query[:30]}...' (Best Score: {best_score:.2f})")
            return default_res

        data = best_match
        pmid = data["Id"] # esummary uses "Id" or we use the list index map. 
        # wait, esummary result is a list of dicts. 'Id' is usually implicitly the Pmid/Id field.
        # Biopython parser: usually has 'Id' key.
        
        # Extract IDs
        doi = ""
        pmcid = ""
        article_ids = data.get("ArticleIds", {})
        
        if "doi" in article_ids:
            doi = article_ids["doi"]
        if "pmc" in article_ids:
            pmcid = article_ids["pmc"]
        if not doi and "DOI" in data:
            doi = data["DOI"]
        if not pmcid and "pmc" in data:
            pmcid = data["pmc"]

        # Formatting
        nlm = format_nlm_citation(data)
        apa = format_apa_citation(data)
        
        # Year
        pub_date = data.get("PubDate", "")
        year = extract_year(pub_date)
        
        # Category
        cat = get_category(data)

        # Link
        link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        
        return {
            "OriginalQuery": query,
            "RefText_NLM": nlm,
            "RefText_APA": apa,
            "PMID": pmid,
            "PMCID": pmcid,
            "DOI": doi,
            "Link": link,
            "Year": year,
            "Category": cat,
            # RIS Fields
            "RIS_Title": data.get("Title", ""),
            "RIS_Journal": data.get("Source", ""),
            "RIS_Volume": data.get("Volume", ""),
            "RIS_Issue": data.get("Issue", ""),
            "RIS_Pages": data.get("Pages", ""),
            "RIS_Authors": data.get("AuthorList", [])
        }

    except Exception as e:
        print(f"Warning: Lookup failed for '{query[:30]}...': {e}")
        return default_res

def clean_fallback_text(text):
    """
    Clean up raw reference text if online lookup fails.
    Removes trailing publisher info like 'LIPPINCOTT WILLIAMS...'.
    Heuristic: Look for standard year;vol:pages pattern and cut after it.
    """
    # Regex for Year;Vol(Issue):Pages (approximate)
    # e.g. 2025;57(10 S):281-281.
    match = re.search(r'(\d{4}\s*;\s*[\w\(\)\s]+:\s*[\w\-]+)', text)
    if match:
        end_idx = match.end()
        # Keep the match, maybe a following dot
        cleaned = text[:end_idx]
        remainder = text[end_idx:]
        if remainder.startswith('.'):
            cleaned += '.'
        return cleaned

    return text

    return text

def parse_fallback_metadata(text):
    """
    Attempt to extract structured data from a raw reference string
    when PubMed lookup fails.
    """
    data = {
        "RIS_Title": "",
        "RIS_Authors": [],
        "RIS_Journal": "",
        "DOI": "",
        "ArXivID": "",
        "Link": ""
    }
    
    # 1. DOI
    doi_match = re.search(r'(10\.\d{4,9}/[-._;()/:A-Za-z0-9]+)', text)
    if doi_match:
        data["DOI"] = doi_match.group(1)
        
    # 2. ArXiv ID
    # Patterns: arXiv:2203.12779 or arXiv preprint arXiv:2506.02260
    arxiv_match = re.search(r'arXiv:(\d+\.\d+)', text, re.IGNORECASE)
    if arxiv_match:
        arxiv_id = arxiv_match.group(1)
        data["ArXivID"] = arxiv_id
        if not data["Link"]:
            data["Link"] = f"https://arxiv.org/abs/{arxiv_id}"
        # If Journal is empty, set it
        data["RIS_Journal"] = "arXiv preprint"

    # 3. Authors and Title (Heuristic)
    # Assumes: "Authors. Title. Journal..."
    parts = text.split('. ')
    if len(parts) >= 2:
        # Part 0 is likely authors
        raw_authors = parts[0]
        # Split authors by comma
        # "Smith J, Doe A" -> ["Smith J", "Doe A"]
        # Filter out things that look like dates or weird stuff? 
        # For now, just simplistic splitting
        auths = [a.strip() for a in raw_authors.split(',')]
        data["RIS_Authors"] = auths
        
        # Part 1 is likely Title
        data["RIS_Title"] = parts[1].strip()
        
        # Part 2 might be Journal, unless it is date
        if len(parts) > 2:
            potential_journal = parts[2].strip()
            # If we didn't set journal yet (e.g. via arXiv detection)
            if not data["RIS_Journal"] and not re.match(r'(19|20)\d{2}', potential_journal):
                data["RIS_Journal"] = potential_journal

    return data

def generate_ris_content(records):
    """
    Generate RIS formatted string from records.
    """
    lines = []
    for rec in records:
        # Start Record
        lines.append("TY  - JOUR")
        
        # Title of main article
        if "RIS_Title" in rec and rec["RIS_Title"]:
             lines.append(f"TI  - {rec['RIS_Title']}")
        else:
             # Fallback title if search failed (heuristic)
             lines.append(f"TI  - {extract_title(rec.get('OriginalQuery', ''))}")

        # Authors
        authors = rec.get("RIS_Authors", [])
        if authors:
            for auth in authors:
                lines.append(f"AU  - {auth}")
        
        # Journal
        if "RIS_Journal" in rec and rec["RIS_Journal"]:
            lines.append(f"JO  - {rec['RIS_Journal']}")
            
        # Year
        if rec["Year"]:
            lines.append(f"PY  - {rec['Year']}")
            
        # Volume
        if "RIS_Volume" in rec and rec["RIS_Volume"]:
            lines.append(f"VL  - {rec['RIS_Volume']}")
            
        # Issue
        if "RIS_Issue" in rec and rec["RIS_Issue"]:
            lines.append(f"IS  - {rec['RIS_Issue']}")
            
        # Pages
        if "RIS_Pages" in rec and rec["RIS_Pages"]:
            # standard can be SP and EP, but generic SP often works
            lines.append(f"SP  - {rec['RIS_Pages']}")
            
        # DOI
        if rec.get("DOI"):
            lines.append(f"DO  - {rec['DOI']}")

        # ArXiv ID (Mapping to AN - Accession Number, and DB - Database Name)
        if rec.get("ArXivID"):
            lines.append(f"AN  - {rec['ArXivID']}")
            lines.append("DB  - arXiv")
            
        # URL
        if rec.get("Link"):
            lines.append(f"UR  - {rec['Link']}")
            
        # PMID (field 'PM' or 'AN' often used)
        if rec.get("PMID"):
             lines.append(f"PM  - {rec['PMID']}") # PubMed ID
        
        lines.append("ER  - ")
        lines.append("") # Empty line between records
        
    return "\n".join(lines)

def build_reference_table(entries):
    """Process entries and fetch info."""
    result = []
    print(f"Processing {len(entries)} references...")
    
    for i, ref in enumerate(entries):
        if (i+1) % 5 == 0 or i == 0:
            print(f"  {i+1}/{len(entries)}...")
            
        title_query = extract_title(ref)
        info = lookup_pubmed_info(title_query, full_ref_text=ref)
        
        # If lookup failed, populate RefText with original
        if not info["PMID"]:
             cleaned_ref = clean_fallback_text(ref)
             info["RefText_NLM"] = cleaned_ref
             info["RefText_APA"] = cleaned_ref
             info["Year"] = extract_year(ref)
             
             # Attempt to parse fallback metadata for RIS
             fallback_data = parse_fallback_metadata(ref)
             info.update(fallback_data)
             
        # Even if lookup SUCCEEDED, we might want to check if we missed DOI/ArXiv in the search results
        # e.g. sometimes PubMed record exists but has no DOI, but the text had it? 
        # For now, let's just prioritize fallback parsing when PMID is missing.
             
        result.append(info)
        
    return result

def save_to_files(records, csv_path, txt_path, style="NLM", sort_by="Newest", group_by_type=True):
    """
    Save to CSV and formatted TXT.
    style: "NLM" or "APA"
    sort_by: "Newest" or "Oldest"
    """
    
    # Sorting
    reverse = (sort_by == "Newest")
    records.sort(key=lambda x: x["Year"], reverse=reverse)
    
    # CSV Write (Combined)
    with open(csv_path, "w", newline='', encoding="utf-8") as f:
        # Added ArXivID to columns
        fieldnames = ["FullCitation", "PMID", "PMCID", "DOI", "ArXivID", "Year", "Link", "Category"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            cit = rec["RefText_NLM"] if style == "NLM" else rec["RefText_APA"]
            writer.writerow({
                "FullCitation": cit,
                "PMID": rec["PMID"],
                "PMCID": rec["PMCID"],
                "DOI": rec["DOI"],
                "ArXivID": rec.get("ArXivID", ""),
                "Link": rec["Link"],
                "Year": rec["Year"],
                "Category": rec["Category"]
            })

    # TXT Write
    with open(txt_path, "w", encoding="utf-8") as f:
        
        # Grouping
        groups = {}
        if group_by_type:
            for rec in records:
                cat = rec["Category"]
                if cat not in groups:
                    groups[cat] = []
                groups[cat].append(rec)
            cat_order = ["Articles", "Abstracts/Posters", "Uncategorized"]
        else:
            groups["All References"] = records
            cat_order = ["All References"]

        first_group = True
        for cat in cat_order:
            if cat not in groups: 
                continue
            
            items = groups[cat]
            if not items: continue
            
            if group_by_type:
                if not first_group:
                    f.write("\n")
                f.write(f"=== {cat} ===\n\n")
            
            for rec in items:
                cit = rec["RefText_NLM"] if style == "NLM" else rec["RefText_APA"]
                parts = [cit]
                
                if rec["DOI"]:
                    parts.append(f"doi: {rec['DOI']}")
                
                if rec.get("ArXivID"):
                    parts.append(f"arXiv: {rec['ArXivID']}")

                if rec["PMID"]:
                    parts.append(f"PMID: {rec['PMID']}")
                
                if rec["PMCID"]:
                    parts.append(f"PMCID: {rec['PMCID']}")
                
                # Add Link at the end
                if rec["Link"]:
                    parts.append(rec["Link"])
                    
                final_str = ". ".join(parts)
                if not final_str.endswith("."):
                    final_str += "."
                final_str = final_str.replace(".. ", ". ").replace("..", ".")
                
                f.write(final_str + "\n\n")
            
            first_group = False
            
    # Also save RIS
    try:
        ris_content = generate_ris_content(records)
        # Use simple name heuristic or passed path. 
        # The CLI call passes csv_path/txt_path. 
        # We'll just replace extension of csv_path for now as a default
        ris_path = csv_path.replace(".csv", ".ris")
        if ris_path == csv_path: ris_path += ".ris"
        
        with open(ris_path, "w", encoding="utf-8") as f:
            f.write(ris_content)
    except Exception as e:
        print(f"Warning: Could not save RIS file: {e}")


if __name__ == "__main__":
    try:
        with open("references.txt", "r", encoding="utf-8") as f:
            raw_text = f.read()
        entries = parse_references(raw_text)
        records = build_reference_table(entries)
        # Default CLI behavior (can act as test)
        save_to_files(records, "recent_references.csv", "recent_references.txt", style="NLM", sort_by="Newest", group_by_type=True)
        print("Done!")
    except FileNotFoundError:
        print("references.txt not found.")
