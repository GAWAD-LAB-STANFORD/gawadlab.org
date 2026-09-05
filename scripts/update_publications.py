#!/usr/bin/env python3
"""Refresh data/publications.json from PubMed, then rebuild the site.

    python3 scripts/update_publications.py

New papers are appended with highlight=false. Existing entries keep any manual
edits (highlight flags, notes). Preprints that later appear as journal articles
must be removed by hand (or add the preprint PMID to SUPERSEDED below).
"""
import json, re, html, time, urllib.request, urllib.parse, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "publications.json"
QUERIES = ['Gawad C[Author]', 'Gawad Charles[Author]']
SUPERSEDED = {"38352480", "39398996"}      # preprints replaced by published versions
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def get(url):
    for attempt in range(6):
        try:
            return urllib.request.urlopen(url, timeout=30).read().decode()
        except Exception as e:
            time.sleep(2 + 2 * attempt)
    sys.exit("PubMed is not responding; try again later.")

def fetch(ids):
    xml = get(BASE + "efetch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids), "retmode": "xml"}))
    out = []
    for art in re.findall(r"<PubmedArticle>.*?</PubmedArticle>", xml, re.S):
        g = lambda pat, default="": (re.search(pat, art, re.S) or [None, default])[1]
        pmid = g(r"<PMID[^>]*>(\d+)")
        title = html.unescape(re.sub(r"<[^>]+>", "", g(r"<ArticleTitle>(.*?)</ArticleTitle>"))).rstrip(".")
        journal = html.unescape(g(r"<Title>(.*?)</Title>"))
        jabbr = html.unescape(g(r"<ISOAbbreviation>(.*?)</ISOAbbreviation>", journal))
        year = g(r"<PubDate>.*?<Year>(\d{4})") or g(r"<MedlineDate>(\d{4})") or "0"
        authors = []
        for a in re.findall(r"<Author [^>]*>(.*?)</Author>", art, re.S):
            ln = re.search(r"<LastName>(.*?)</LastName>", a); ini = re.search(r"<Initials>(.*?)</Initials>", a)
            if ln: authors.append(html.unescape(ln.group(1)) + " " + (ini.group(1) if ini else ""))
        out.append({"pmid": pmid, "doi": g(r'<ArticleId IdType="doi">(.*?)</ArticleId>'), "pmc": g(r'<ArticleId IdType="pmc">(.*?)</ArticleId>'),
                    "year": int(year), "title": title, "journal": jabbr, "journal_full": journal, "authors": authors,
                    "preprint": jabbr in ("bioRxiv", "medRxiv"), "highlight": False, "note": ""})
    return out

def main():
    existing = json.load(open(DATA)) if DATA.exists() else []
    known = {p["pmid"] for p in existing}
    ids = set()
    for q in QUERIES:
        r = json.loads(get(BASE + "esearch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "term": q, "retmax": 500, "retmode": "json"})))
        ids |= set(r["esearchresult"]["idlist"]); time.sleep(0.5)
    new_ids = sorted(ids - known - SUPERSEDED, key=int)
    if new_ids:
        added = fetch(new_ids)
        for p in added: print("NEW:", p["year"], p["journal"], "-", p["title"][:80])
        existing += added
    else:
        print("No new PubMed records.")
    existing.sort(key=lambda p: (-p["year"], -int(p["pmid"])))
    json.dump(existing, open(DATA, "w"), indent=1, ensure_ascii=False)
    print(f"{len(existing)} publications in {DATA.relative_to(ROOT)}")
    print("Check new entries for other authors named 'C Gawad' before publishing, then set highlight=true on any you want featured.")
    subprocess.run([sys.executable, str(ROOT / "build.py")], check=True)

if __name__ == "__main__":
    main()
