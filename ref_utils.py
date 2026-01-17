import re
import csv
import time
import ssl
from datetime import datetime

print("RefCheck module v2 loaded")

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


def parse_references_text(text):
    """
    Split the input text (TXT) into references.
    """
    # Return dict structure
    entries = [{"title": None, "text": line.strip(), "metadata": {}} for line in text.splitlines() if line.strip()]
    return entries

def parse_csv_content(text):
    """
    Parse CSV content string. Expects 'Title' column.
    """
    import io
    entries = []
    try:
        reader = csv.DictReader(io.StringIO(text))
        for row in reader:
            # Flexible column checking
            if "Title" in row and row["Title"]:
                title = row["Title"].strip()
                # Handle potential BOM in Authors column
                authors = row.get("Authors") or row.get("Author") or row.get("\ufeffAuthors") or ""
                
                # Extract other metadata for fallback
                year = row.get("Year", "").strip()
                source = row.get("Publication") or row.get("Journal") or row.get("Booktitle") or ""
                volume = row.get("Volume", "").strip()
                issue = row.get("Number", "").strip()
                pages = row.get("Pages", "").strip()
                
                # Construct full text for display/fallback context
                # We try to mimic a standard citation string for the 'text' field
                cit_details = year
                if volume: cit_details += f";{volume}"
                if issue: cit_details += f"({issue})"
                if pages: cit_details += f":{pages}"

                full_text_parts = [authors, title, source, cit_details]
                full_text = ". ".join([p for p in full_text_parts if p]) + "."
                
                entries.append({
                    "title": title,
                    "text": full_text,
                    "metadata": {
                        "Year": year,
                        "Source": source,
                        "Volume": volume,
                        "Issue": issue,
                        "Pages": pages,
                        "Authors": authors
                    }
                })
            elif "Citation" in row and row["Citation"]:
                 entries.append({"title": None, "text": row["Citation"], "metadata": {}})
    except Exception as e:
        print(f"Error parsing CSV: {e}")
    return entries

def format_bibtex_authors(bibtex_author_str):
    """
    Format BibTeX 'Last, First and Last, First' to 'Last F, Last F'.
    Also handles 'others' -> 'et al'.
    """
    if not bibtex_author_str: return ""
    
    # Split by ' and ' (BibTeX delimiter)
    # Check for ' and ' - some might be single author
    authors = bibtex_author_str.split(' and ')
    formatted = []
    
    for auth in authors:
        clean_auth = auth.strip()
        if clean_auth.lower() == 'others':
            formatted.append("et al")
            continue
            
        # Parse "Last, First"
        if ',' in clean_auth:
            parts = clean_auth.split(',')
            last = parts[0].strip()
            first = parts[1].strip() if len(parts) > 1 else ""
            
            # Initials only? Or keep full first?
            # User wants "Display" like PubMed: "Last FM"
            # But let's do "Last F" for simplicity and robustness
            initials = "".join([n[0] for n in first.split() if n])
            if initials:
                formatted.append(f"{last} {initials}")
            else:
                formatted.append(last)
        else:
            # Maybe "First Last"? Assume string is name.
            formatted.append(clean_auth)
            
    return ", ".join(formatted)

def clean_latex(text):
    """
    Remove LaTeX formatting for plain text use.
    """
    if not text: return ""
    # Unescape: \& -> &, \_ -> _
    text = text.replace(r"\&", "&").replace(r"\_", "_").replace(r"\%", "%")
    # Dashes: -- -> -
    text = text.replace("--", "-")
    # Braces: {Title} -> Title. Be simple: remove all braces.
    text = text.replace("{", "").replace("}", "")
    # Extra whitespace cleanup
    return " ".join(text.split())

def parse_bibtex_content(text):
    """
    Parse BibTeX content string using Regex (Simple).
    """
    entries = []
    
    # Heuristic: split by @, then look for key fields
    raw_entries = text.split('@')
    for raw in raw_entries:
        if not raw.strip(): continue
        
        # Helper to extract field
        def get_field(name):
            match = re.search(rf'{name}\s*=\s*[\"{{](.+?)[\"}}]\s*[,}}]', raw, re.IGNORECASE | re.DOTALL)
            val = " ".join(match.group(1).split()) if match else ""
            return clean_latex(val)

        title = get_field("title")
        if title:
            # Note: Author cleaning is handled partly by format_bibtex_authors but good to have clean input
            # However format_bibtex_authors expects raw structure ' and ' so be careful not to break it
            # But ' and ' is standard English, LaTeX chars inside names?
            # Safe to use clean_latex on author string provided ' and ' isn't ' \& ' (unlikely)
            raw_authors = get_field("author")
            
            # Use raw uncleaned for formatting logic?
            # Actually get_field now cleans. 
            # If BibTeX had "Smith \& Wesson", clean_latex makes it "Smith & Wesson".
            # format_bibtex_authors splits by " and ". 
            # If author was "Barnes \and Noble", regex matching might be tricky but usually " and ".
            # Let's assume standard ' and '.
            
            # Format nicely for display
            authors = format_bibtex_authors(raw_authors)
            
            year = get_field("year")
            # checked for journal OR booktitle
            source = get_field("journal") or get_field("booktitle")
            volume = get_field("volume")
            number = get_field("number")
            pages = get_field("pages")
            
            # Construct full text for display/validation
            # Format: 'Authors. Title. Journal. Year;Volume(Issue):Pages.'
            cit_details = year
            if volume: cit_details += f";{volume}"
            if number: cit_details += f"({number})"
            if pages: cit_details += f":{pages}"

            full_text_parts = [authors, title, source, cit_details]
            full_text = ". ".join([p for p in full_text_parts if p]) + "."
            
            entries.append({
                "title": title,
                "text": full_text,
                "metadata": {
                    "Year": year,
                    "Source": source,
                    "Volume": volume,
                    "Issue": number,
                    "Pages": pages,
                    "Authors": authors
                }
            })
                
    return entries

def parse_any_input(content, filename):
    """
    Dispatch based on filename extension.
    """
    ext = filename.lower().split('.')[-1]
    if ext == 'csv':
        return parse_csv_content(content)
    elif ext in ['bib', 'bibtex']:
        return parse_bibtex_content(content)
    else:
        return parse_references_text(content)

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
    Determine category: Article, Abstract, Preprint, or Other.
    """
    # Check for Preprints first
    source = item.get("Source", "").lower()
    if "arxiv" in source or "medrxiv" in source or "biorxiv" in source:
        return "Preprint"
        
    ptypes = item.get("PubTypeList", [])
    ptypes_lower = [p.lower() for p in ptypes]
    if any("abstract" in p for p in ptypes_lower):
        return "Abstracts/Posters"
    return "Articles" 

def lookup_pubmed_info(query, full_ref_text=None, email="A.N.Other@example.com"):
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

    Entrez.email = email
    
    try:
        # Build Composite Query
        qs = [f"({query})"]
        
        # 1. Extracted Title Query
        # Try to extract a clean title from the query/ref_text
        text_to_use = full_ref_text if full_ref_text else query
        extracted_title = extract_title(text_to_use)
        
        if extracted_title and len(extracted_title) > 10:
             qs.append(f"({extracted_title}[Title])")
             
        # 2. Smart Keyword Query (Auth + Title Keywords)
        if full_ref_text:
            try:
                # 1. BibTeX style " and "
                if " and " in full_ref_text[:200]: 
                    auths = full_ref_text.split('.')
                    if len(auths) > 0:
                        raw_auths = auths[0].split(' and ')
                        surnames = [a.strip().split(',')[0].strip().split(' ')[0] for a in raw_auths[:4]]
                else: 
                    # 2. Comma style
                    potential_list = full_ref_text.split('.')[0]
                    raw_parts = potential_list.split(',')
                    surnames = []
                    for p in raw_parts[:4]:
                        s = p.strip().split(' ')[0]
                        if len(s) > 1 and not s.isdigit():
                            surnames.append(s)
            except:
                surnames = []

            if extracted_title and len(extracted_title) > 10 and surnames:
                keywords = " ".join(extracted_title.split()[:6])
                auth_query = " OR ".join([f"{s}[Author]" for s in surnames if len(s)>1])
                if auth_query:
                    qs.append(f"(({auth_query}) AND {keywords}[Title])")
        
        # Join into single composite query
        composite_query = " OR ".join(qs)
        
        # 1. Single API Call
        # Removed sleep for performance as per user request
        search_handle = Entrez.esearch(db="pubmed", term=composite_query, retmax=20) 
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_results["IdList"]

        if not id_list:
            return default_res
            
        # 2. Fetch Summaries for candidates
        # Removed sleep
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
            # Updated: >90% correlation required
            similarity = SequenceMatcher(None, clean_found, clean_query).ratio()
            
            # Also check substring match (sometimes query is cut off or data has extra info)
            if clean_query in clean_found or clean_found in clean_query:
                # Boost if substring match exists
                if similarity < 0.92: 
                    similarity = 0.92 # Boost to pass strict check if substring matches
            
            if similarity < 0.90: # Very strict threshold
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
        
        # Check if authors are separated by " and " (BibTeX style) or "," (Standard)
        if " and " in raw_authors:
             # Split by " and ", and cleanup
             auths = [a.strip() for a in raw_authors.split(' and ')]
        else:
             # Standard comma split
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

    # 4. Citation Details (Year, Volume, Issue, Pages)
    # Pattern: 2024;61(4):1876-1887 or 2024. 61(4):1876-1887
    # Updated: make pages optional
    cit_match = re.search(r'(\d{4})[;.]\s*(\d+)[\(](\d+)[\)](?::\s*(\d+(?:-\d+)?))?', text)
    if cit_match:
        # Group 1: Year, 2: Vol, 3: Issue, 4: Pages (optional)
        data["RIS_Volume"] = cit_match.group(2)
        data["RIS_Issue"] = cit_match.group(3)
        if cit_match.group(4):
            data["RIS_Pages"] = cit_match.group(4)
    else:
        # Try generic pattern Year;Vol:Pages (no issue)
        cit_match_2 = re.search(r'(\d{4})[;.]\s*(\d+):\s*(\d+(?:-\d+)?)', text)
        if cit_match_2:
             data["RIS_Volume"] = cit_match_2.group(2)
             data["RIS_Pages"] = cit_match_2.group(3)

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



def save_to_files(records, csv_path, txt_path, style="NLM", sort_by="Newest", group_by_type=True):
    """
    Save to CSV and formatted TXT.
    style: "NLM" or "APA"
    sort_by: "Newest" or "Oldest"
    """
    
    # Sorting
    reverse = (sort_by == "Newest")
    # robust sort key
    records.sort(key=lambda x: int(x["Year"]) if str(x["Year"]).isdigit() else 0, reverse=reverse)
    
    # CSV Write (Combined)
    with open(csv_path, "w", newline='', encoding="utf-8-sig") as f:
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
    import argparse
    import os
    
    parser = argparse.ArgumentParser(description="Process reference files.")
    parser.add_argument("input_dir", nargs='?', default="input_example", help="Input directory containing .csv, .txt, or .bib files")
    parser.add_argument("output_dir", nargs='?', default="output", help="Output directory for results")
    
    args = parser.parse_args()
    
    input_dir = args.input_dir
    output_dir = args.output_dir
    
    try:        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        supported_exts = ['.txt', '.csv', '.bib', '.bibtex']
        
        if not os.path.exists(input_dir):
            print(f"Input directory '{input_dir}' not found.")
        else:
            files_found = [f for f in os.listdir(input_dir) if any(f.endswith(ext) for ext in supported_exts)]
            # Filter specifically for what user might want if many files exist? 
            # No, just process all supported files in dir.
            print(f"Found {len(files_found)} files to process in {input_dir}: {files_found}")
            
            for filename in files_found:
                file_path = os.path.join(input_dir, filename)
                print(f"\n--- Processing {filename} ---")
                
                with open(file_path, "r", encoding="utf-8", errors='replace') as f:
                    content = f.read()
                
                entries = parse_any_input(content, filename)
                if not entries:
                    print(f"  No entries found in {filename}.")
                    continue
                    
                records = []
                print(f"  Found {len(entries)} references/entries...")
                for i, entry in enumerate(entries):
                    if (i+1) % 5 == 0 or i == 0:
                         print(f"  {i+1}/{len(entries)}...")
                    
                    # Entry is a dict: {'title': ..., 'text': ...}
                    ref_text = entry['text']
                    
                    # Use explicit title if available, otherwise extract
                    if entry['title']:
                        t = entry['title']
                    else:
                        t = extract_title(ref_text)
                        
                    # Pass the full ref string (which might be author+title) for validation
                    info = lookup_pubmed_info(t, full_ref_text=ref_text)
                    
                    if not info["PMID"]:
                        cleaned = clean_fallback_text(ref_text)
                        info["RefText_NLM"] = cleaned
                        info["RefText_APA"] = cleaned
                        
                        # Use metadata from parser if available (CSV/BibTeX), else extract from text
                        meta = entry.get('metadata', {})
                        info["Year"] = meta.get('Year') or extract_year(ref_text)
                        
                        # Update fallback data with parser metadata
                        # We merge: parser metadata > extracted metadata (via fallback parser)
                        fallback_extracted = parse_fallback_metadata(ref_text)
                        info.update(fallback_extracted)
                        
                        if meta.get('Source'): 
                            info['Journal'] = meta['Source']
                            info['RIS_Journal'] = meta['Source']
                        if meta.get('Volume'): 
                            info['Volume'] = meta['Volume']
                            info['RIS_Volume'] = meta['Volume']
                        if meta.get('Issue'): 
                            info['Issue'] = meta['Issue']
                            info['RIS_Issue'] = meta['Issue']
                        if meta.get('Pages'): 
                            info['Pages'] = meta['Pages']
                            info['RIS_Pages'] = meta['Pages']
                            info['RIS_Issue'] = meta['Issue']
                        if meta.get('Pages'): 
                            info['Pages'] = meta['Pages']
                            info['RIS_Pages'] = meta['Pages']
                        
                        # Populate RIS Authors from metadata if available to avoid bad re-parsing
                        if meta.get('Authors'):
                            raw_auths = meta['Authors']
                            if " and " in raw_auths:
                                info['RIS_Authors'] = [a.strip() for a in raw_auths.split(' and ')]
                            else:
                                # Fallback or already comma separated?
                                # If CSV, likely comma or semicolon. Let's assume comma for now or no split if complex.
                                info['RIS_Authors'] = [a.strip() for a in raw_auths.split(',')]
                        
                        # Re-format text with new metadata
                        # (Simple re-construction if we have good metadata)
                        if info['Journal'] and info['Year']:
                             # NLM Style: Authors. Title. Journal. Year;Vol(Issue):Pages.
                             # This is rough, but better than nothing.
                             pass
                        
                    records.append(info)
                
                # Derive output name
                base_name = os.path.splitext(filename)[0]
                out_csv = os.path.join(output_dir, f"{base_name}_results.csv")
                out_txt = os.path.join(output_dir, f"{base_name}_results.txt")
                
                save_to_files(records, out_csv, out_txt, style="NLM", sort_by="Newest", group_by_type=True)
                print(f"  Saved results to {out_csv}")

        print("\nBatch processing complete!")
    except Exception as e:
        print(f"An error occurred: {e}")
