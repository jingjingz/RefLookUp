import streamlit as st
import pandas as pd
import io
import csv
import refcheck

st.set_page_config(page_title="RefLookUp", page_icon="📚", layout="wide")

st.title("📚 RefLookUp: Reference Metadata Search")
st.markdown("""
Upload a **text file** containing your references (one per line). 
The app will search PubMed for PMIDs, PMCIDs, DOIs, and updated citation metadata.
""")

# Initialize Session State
if "results" not in st.session_state:
    st.session_state["results"] = None

# Sidebar Options
st.sidebar.header("Options")
# Email for Entrez
entrez_email = st.sidebar.text_input("Entrez Email (Required for PubMed)", value="", placeholder="your.email@example.com", help="NCBI requires a valid email for API access. Providing this prevents rate-limiting errors.")

style_option = st.sidebar.selectbox("Citation Style", ["NLM", "APA"])
sort_option = st.sidebar.radio("Sort Order", ["Newest", "Oldest"])
group_option = st.sidebar.checkbox("Group by Type (Articles/Abstracts)", value=True)

uploaded_file = st.file_uploader("Choose a file (references.txt)", type=["txt"])

if uploaded_file is not None:
    # Button to trigger search
    if st.button("Search References"):
        if not entrez_email or "@" not in entrez_email:
             st.warning("⚠️ Please enter a valid email address in the sidebar to ensure PubMed access.")
        
        stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
        raw_text = stringio.read()
        
        # Use provided email or fallback (though fallback often fails on cloud)
        email_to_use = entrez_email if entrez_email else "A.N.Other@example.com"
        
        with st.spinner("Searching PubMed... This may take a moment."):
            entries = refcheck.parse_references(raw_text)
            results = refcheck.build_reference_table(entries, email=email_to_use)
            st.session_state["results"] = results
            st.success(f"Processed {len(results)} references.")

# Display Results if they exist in state
if st.session_state["results"]:
    results = st.session_state["results"]
    
    # Sorting
    reverse_sort = (sort_option == "Newest")
    # Sort a copy for display/export
    sorted_results = sorted(results, key=lambda x: x["Year"], reverse=reverse_sort)
    
    # Prepare DataFrame
    df = pd.DataFrame(sorted_results)
    
    # Ensure Citation Column based on style
    cit_col = "RefText_NLM" if style_option == "NLM" else "RefText_APA"
    df["FullCitation"] = df[cit_col]
    
    # Display Columns
    # Added ArXivID
    disp_cols = ["FullCitation", "PMID", "PMCID", "DOI", "ArXivID", "Year", "Link", "Category"]
    
    # Ensure columns exist in DF
    for c in disp_cols:
        if c not in df.columns: df[c] = ""
        
    st.subheader("Results Table")
    
    # Config for clickable links
    column_config = {
        "Link": st.column_config.LinkColumn(
            "Link",
            help="Click to open reference",
            validate="^https://.*",
            display_text="View Paper"
        )
    }
    
    st.dataframe(
        df[disp_cols],
        width="stretch",
        column_config=column_config,
        hide_index=True
    )
    
    st.divider()
    st.subheader("Download Results")
    
    c1, c2, c3 = st.columns(3)
    
    # 1. CSV Download
    csv_buffer = io.StringIO()
    writer = csv.DictWriter(csv_buffer, fieldnames=disp_cols)
    writer.writeheader()
    for rec in sorted_results:
        cit = rec["RefText_NLM"] if style_option == "NLM" else rec["RefText_APA"]
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
        
    with c1:
        st.download_button(
            label="Download CSV",
            data=csv_buffer.getvalue(),
            file_name="recent_references.csv",
            mime="text/csv"
        )
    
    # 2. TXT Download
    txt_buffer = io.StringIO()
    
    # Grouping Logic
    groups = {}
    if group_option:
        for r in sorted_results:
            cat = r["Category"]
            if cat not in groups: groups[cat] = []
            groups[cat].append(r)
        cat_order = ["Articles", "Abstracts/Posters", "Uncategorized"]
    else:
        groups["All References"] = sorted_results
        cat_order = ["All References"]
    
    first_group = True
    for cat in cat_order:
        if cat not in groups: continue
        items = groups[cat]
        if not items: continue
        
        if group_option:
            if not first_group: txt_buffer.write("\n")
            txt_buffer.write(f"=== {cat} ===\n\n")
        
        for rec in items:
            cit = rec["RefText_NLM"] if style_option == "NLM" else rec["RefText_APA"]
            parts = [cit]
            
            if rec["DOI"]: parts.append(f"doi: {rec['DOI']}")
            if rec.get("ArXivID"): parts.append(f"arXiv: {rec['ArXivID']}")
            if rec["PMID"]: parts.append(f"PMID: {rec['PMID']}")
            if rec["PMCID"]: parts.append(f"PMCID: {rec['PMCID']}")
            
            # Link at the end
            if rec.get("Link"): parts.append(rec["Link"])
            
            final_str = ". ".join(parts)
            # Remove double dots if any
            final_str = final_str.replace(".. ", ". ").replace("..", ".")
            # Ensure it ends with dot? URLs shouldn't strictly end with dot if it breaks them, 
            # but standard citation style often does.
            # If the last part is a Link, verify avoiding dot if preferred, 
            # but user complaint was 'missing', not 'broken'. 
            # We'll stick to 'ends with dot' for consistency unless it's a raw URL line.
            if not final_str.endswith("."): final_str += "."
            
            txt_buffer.write(final_str + "\n\n")
        
        first_group = False
    
    with c2:
        st.download_button(
            label="Download Text File",
            data=txt_buffer.getvalue(),
            file_name="recent_references.txt",
            mime="text/plain"
        )

    # 3. RIS Download
    ris_content = refcheck.generate_ris_content(sorted_results)
    with c3:
        st.download_button(
            label="Download RIS (NCBI)",
            data=ris_content,
            file_name="recent_references.ris",
            mime="application/x-research-info-systems"
        )
