#!/usr/bin/env python3
"""Build the Gawad Lab static site.

    python3 build.py            # writes the finished site into ./site/

Pages are HTML fragments in src/pages/. The first line of each fragment is a
comment with the page title and description. Placeholders such as
{{LATEST_PUBS}} are filled from the JSON files in data/. publications.html is
generated entirely from data/publications.json.
"""
import html, json, re, shutil, datetime, hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC, DATA, OUT = ROOT / "src", ROOT / "data", ROOT / "site"
SITE_URL = "https://gawadlab.org"
YEAR = datetime.date.today().year
def _v(rel):
    """Short content hash so browsers fetch changed CSS/JS immediately."""
    return hashlib.md5((SRC / rel).read_bytes()).hexdigest()[:8]

# Lab members whose names get bolded in author lists (surname + initials as PubMed writes them).
LAB_AUTHORS = {"Gawad C", "Gonzalez-Pena V", "Pang Y", "Wardhani K", "Klein D", "Aragon A", "Cai B",
               "Schulz S", "Xia Y", "Natarajan S", "Carter RA", "Carter R", "Mahmud O", "Inaba Y",
               "Youssef S", "Agarwal V", "Cherian A"}

NAV = [("index.html", "Home"), ("research.html", "Research"), ("people.html", "People"),
       ("publications.html", "Publications"), ("news.html", "News"), ("resources.html", "Resources")]

ICONS = {
    "mail": '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><rect x="3" y="5" width="18" height="14" rx="2"/><path d="m3 7 9 6 9-6"/></svg>',
    "github": '<svg viewBox="0 0 24 24" fill="currentColor"><path d="M12 .5C5.7.5.5 5.7.5 12c0 5.1 3.3 9.4 7.9 10.9.6.1.8-.3.8-.6v-2c-3.2.7-3.9-1.4-3.9-1.4-.5-1.3-1.3-1.7-1.3-1.7-1-.7.1-.7.1-.7 1.2.1 1.8 1.2 1.8 1.2 1 1.8 2.7 1.3 3.4 1 .1-.8.4-1.3.7-1.6-2.6-.3-5.3-1.3-5.3-5.7 0-1.3.5-2.3 1.2-3.1-.1-.3-.5-1.5.1-3.1 0 0 1-.3 3.2 1.2a11 11 0 0 1 5.8 0c2.2-1.5 3.2-1.2 3.2-1.2.6 1.6.2 2.8.1 3.1.8.8 1.2 1.8 1.2 3.1 0 4.4-2.7 5.4-5.3 5.7.4.4.8 1.1.8 2.2v3.2c0 .3.2.7.8.6 4.6-1.5 7.9-5.8 7.9-10.9C23.5 5.7 18.3.5 12 .5z"/></svg>',
    "scholar": '<svg viewBox="0 0 24 24" fill="currentColor"><path d="M12 3 1 9l11 6 9-4.9V17h2V9L12 3zm-6 10.6V17c0 1.7 2.7 3.5 6 3.5s6-1.8 6-3.5v-3.4l-6 3.3-6-3.3z"/></svg>',
    "linkedin": '<svg viewBox="0 0 24 24" fill="currentColor"><path d="M4.98 3.5C4.98 4.9 3.9 6 2.5 6S0 4.9 0 3.5 1.1 1 2.5 1s2.48 1.1 2.48 2.5zM.2 8h4.6v14H.2V8zm7.6 0h4.4v1.9h.1c.6-1.1 2.1-2.3 4.3-2.3 4.6 0 5.5 3 5.5 7V22h-4.6v-6.6c0-1.6 0-3.6-2.2-3.6s-2.5 1.7-2.5 3.5V22H7.8V8z"/></svg>',
}

def esc(s): return html.escape(str(s), quote=True)

def head(title, desc, page, og_image="assets/og-image.jpg"):
    full = f"{title} · Gawad Lab" if page != "index.html" else "Gawad Lab · Stanford Medicine"
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(full)}</title>
<meta name="description" content="{esc(desc)}">
<link rel="canonical" href="{SITE_URL}/{'' if page=='index.html' else page}">
<meta property="og:type" content="website">
<meta property="og:site_name" content="Gawad Lab">
<meta property="og:title" content="{esc(full)}">
<meta property="og:description" content="{esc(desc)}">
<meta property="og:url" content="{SITE_URL}/{'' if page=='index.html' else page}">
<meta property="og:image" content="{SITE_URL}/{og_image}">
<meta name="twitter:card" content="summary_large_image">
<link rel="icon" href="assets/favicon.svg" type="image/svg+xml">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:ital,opsz,wght,SOFT,WONK@0,9..144,400..700,0..100,0..1;1,9..144,400..700,0..100,0..1&family=Inter:wght@400;500;600;700;800&display=swap" rel="stylesheet">
<link rel="stylesheet" href="css/style.css?v={_v("css/style.css")}">
<script type="application/ld+json">{{"@context":"https://schema.org","@type":"ResearchOrganization","name":"Gawad Lab","url":"{SITE_URL}","parentOrganization":{{"@type":"CollegeOrUniversity","name":"Stanford University School of Medicine"}},"address":{{"@type":"PostalAddress","streetAddress":"240 Pasteur Dr, BMI Building RM2200","addressLocality":"Palo Alto","addressRegion":"CA","postalCode":"94304","addressCountry":"US"}},"email":"cgawad@stanford.edu"}}</script>
</head>
<body>
<a class="skip-link" href="#main">Skip to content</a>
<header class="site-header">
  <div class="container nav">
    <a class="brand" href="index.html" aria-label="Gawad Lab home"><img src="assets/logo.png" alt="Stanford Medicine · Gawad Lab" width="592" height="164"></a>
    <button class="nav-toggle" aria-expanded="false" aria-controls="nav-links" aria-label="Menu"><span></span></button>
    <ul class="nav-links" id="nav-links">
""" + "".join(f'      <li><a href="{h}"{" class=\"active\"" if h == page else ""}>{t}</a></li>\n' for h, t in NAV) + """      <li class="nav-cta"><a href="join.html">Join the lab</a></li>
    </ul>
  </div>
</header>
<main id="main">
"""

def footer():
    return f"""</main>
<footer class="site-footer">
  <div class="container">
    <div class="footer-grid">
      <div class="footer-brand">
        <img src="assets/logo.png" alt="Stanford Medicine · Gawad Lab">
        <p>We invent single-cell and cell-free genomics technologies and use them to understand how childhood cancers arise, evolve, and resist treatment.</p>
        <div class="social">
          <a href="https://github.com/GAWAD-LAB-STANFORD" aria-label="GitHub">{ICONS['github']}</a>
          <a href="https://scholar.google.com/citations?hl=en&user=Nbk0c_oAAAAJ" aria-label="Google Scholar">{ICONS['scholar']}</a>
          <a href="mailto:cgawad@stanford.edu" aria-label="Email">{ICONS['mail']}</a>
        </div>
      </div>
      <div>
        <h4>Explore</h4>
        <ul>{"".join(f'<li><a href="{h}">{t}</a></li>' for h, t in NAV)}<li><a href="join.html">Join the lab</a></li></ul>
      </div>
      <div>
        <h4>Contact</h4>
        <address>Charles Gawad, MD, PhD<br>Department of Pediatrics, Stanford Medicine<br>240 Pasteur Dr, BMI Building RM2200<br>Palo Alto, CA 94304<br><a href="mailto:cgawad@stanford.edu">cgawad@stanford.edu</a></address>
      </div>
    </div>
    <div class="footer-bottom">
      <span>© {YEAR} Gawad Lab, Stanford University.</span>
      <span>Cover art from the lab's PNAS 2021 issue · original artwork by SciStories.</span>
    </div>
  </div>
</footer>
<script src="js/main.js?v={_v("js/main.js")}"></script>
</body>
</html>
"""

# ---------- publications ----------
def author_html(authors):
    def fmt(a):
        return f"<b>{esc(a)}</b>" if a in LAB_AUTHORS else esc(a)
    n = len(authors)
    if n <= 10:
        return ", ".join(fmt(a) for a in authors)
    head_, tail = authors[:6], authors[6:]
    return (", ".join(fmt(a) for a in head_) + '<span class="ellipsis">, … ' + fmt(authors[-1]) + "</span>"
            + '<span class="more">, ' + ", ".join(fmt(a) for a in tail) + "</span> "
            + f'<button class="link-btn" type="button" data-count="{n}">show all {n} authors</button>')

def pub_item(p, compact=False):
    doi = f"https://doi.org/{p['doi']}" if p.get("doi") else f"https://pubmed.ncbi.nlm.nih.gov/{p['pmid']}/"
    lab = "1" if any(a in LAB_AUTHORS for a in p["authors"][:1] + p["authors"][-1:]) else "0"
    links = [f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p["pmid"]}/">PubMed</a>']
    if p.get("doi"): links.append(f'<a href="https://doi.org/{esc(p["doi"])}">DOI</a>')
    if p.get("pmc"): links.append(f'<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/{p["pmc"]}/">Free full text</a>')
    badge = ' <span class="badge">Preprint</span>' if p.get("preprint") else ""
    search = " ".join([p["title"], p["journal"], str(p["year"])] + p["authors"]).lower()
    cls = "pub" + (" highlight" if p.get("highlight") else "")
    authors = author_html(p["authors"]) if not compact else ", ".join(esc(a) for a in p["authors"][:3]) + (" et al." if len(p["authors"]) > 3 else "")
    return f"""<li class="{cls}" data-search="{esc(search)}" data-year="{p['year']}" data-preprint="{'1' if p.get('preprint') else '0'}" data-lab="{lab}">
  <h3 class="pub-title"><a href="{esc(doi)}">{esc(p['title'])}</a></h3>
  <p class="pub-authors">{authors}</p>
  <div class="pub-meta"><span class="journal">{esc(p['journal'])}</span><span>{p['year']}</span>{badge}</div>
  {"" if compact else '<div class="pub-links">' + " ".join(links) + "</div>"}
</li>"""

def publications_page(pubs):
    body = [f"""<section class="page-hero"><div class="container">
  <span class="eyebrow">Publications</span>
  <h1>Papers from the lab and with our <em class="accent">collaborators</em></h1>
  <p class="lede">Peer-reviewed articles and preprints, newest first. Lab members are shown in bold. The list is generated from PubMed; see <a href="https://scholar.google.com/citations?hl=en&user=Nbk0c_oAAAAJ">Google Scholar</a> for citation metrics.</p>
</div></section>
<section class="section"><div class="container">
  <div class="pub-tools">
    <label class="visually-hidden" for="pub-search">Search publications</label>
    <input id="pub-search" type="search" placeholder="Search by title, author, journal, or year…" autocomplete="off">
    <span class="count" id="pub-count"></span>
  </div>
  <div class="chip-row" style="margin-bottom:1.5rem">
    <button class="chip" data-filter="all" aria-pressed="true">All</button>
    <button class="chip" data-filter="highlight" aria-pressed="false">Selected</button>
    <button class="chip" data-filter="lab" aria-pressed="false">Lab-led</button>
    <button class="chip" data-filter="preprint" aria-pressed="false">Preprints</button>
  </div>
"""]
    current = None
    for p in pubs:
        if p["year"] != current:
            if current is not None: body.append("</ul>")
            current = p["year"]
            body.append(f'<h2 class="pub-year" data-year="{current}">{current}</h2>\n<ul class="pub-list">')
        body.append(pub_item(p))
    body.append('</ul>\n<p class="no-results" id="pub-empty" hidden>No publications match your search.</p>\n</div></section>')
    return "\n".join(body)

# ---------- people ----------
def person_card(m):
    links = ""
    if m.get("links"):
        L = m["links"]; parts = []
        if L.get("email"): parts.append(f'<a href="mailto:{esc(L["email"])}" aria-label="Email">{ICONS["mail"]}</a>')
        if L.get("scholar"): parts.append(f'<a href="{esc(L["scholar"])}" aria-label="Google Scholar">{ICONS["scholar"]}</a>')
        if L.get("github"): parts.append(f'<a href="{esc(L["github"])}" aria-label="GitHub">{ICONS["github"]}</a>')
        if L.get("linkedin"): parts.append(f'<a href="{esc(L["linkedin"])}" aria-label="LinkedIn">{ICONS["linkedin"]}</a>')
        links = '<div class="links">' + "".join(parts) + "</div>"
    return f"""<article class="person">
  <div class="portrait"><img src="assets/people/{esc(m['image'])}" alt="{esc(m['name'])}" width="150" height="150" loading="lazy"></div>
  <h3>{esc(m['name'])}</h3>
  <p class="role">{esc(m['role'])}</p>
  <p class="bio">{esc(m['bio'])}</p>{links}
</article>"""

def alum_card(m):
    yrs = f" · {esc(m['years'])}" if m.get("years") else ""
    return f"""<article class="alum">
  <img src="assets/people/{esc(m['image'])}" alt="" width="64" height="64" loading="lazy">
  <div><h3>{esc(m['name'])}</h3><p class="role">{esc(m['role'])}{yrs}</p><p class="bio">{esc(m['bio'])}</p></div>
</article>"""

def news_items(items):
    out = []
    for it in items:
        title = f'<a href="{esc(it["link"])}">{esc(it["title"])}</a>' if it.get("link") else esc(it["title"])
        out.append(f'<li><time datetime="{esc(it["date"])}">{esc(it["date"])}</time><h3>{title}</h3><p>{esc(it["text"])}</p></li>')
    return "\n".join(out)

def resource_groups(groups):
    out = []
    for g in groups:
        lis = "".join(f'<li><a href="{esc(u)}" rel="noopener">{esc(t)}</a></li>' for t, u in g["links"])
        out.append(f'<div class="resource-group"><h3>{esc(g["title"])}</h3><ul>{lis}</ul></div>')
    return "\n".join(out)

FUNDERS = [("alsf.jpg", "Alex's Lemonade Stand Foundation"), ("bwf.png", "Burroughs Wellcome Fund"),
           ("hyundai-hope-on-wheels.jpg", "Hyundai Hope On Wheels"), ("nih-new-innovator.png", "NIH Director's New Innovator Award"),
           ("stanford.png", "Stanford University"), ("ash.png", "American Society of Hematology"),
           ("cz-biohub.png", "Chan Zuckerberg Biohub"), ("lls.png", "Leukemia & Lymphoma Society")]

def build():
    pubs = json.load(open(DATA / "publications.json"))
    people = json.load(open(DATA / "people.json"))
    news = json.load(open(DATA / "news.json"))["items"]
    resources = json.load(open(DATA / "resources.json"))["groups"]

    if OUT.exists(): shutil.rmtree(OUT)
    OUT.mkdir()
    for d in ("css", "js", "assets"): shutil.copytree(SRC / d, OUT / d)

    fills = {
        "{{LATEST_PUBS}}": "\n".join(pub_item(p, compact=True) for p in pubs[:4]),
        "{{PEOPLE_CURRENT}}": "\n".join(person_card(m) for m in people["current"]),
        "{{PEOPLE_ALUMNI}}": "\n".join(alum_card(m) for m in people["alumni"]),
        "{{PEOPLE_MASCOTS}}": "\n".join(f'<div class="mascot"><img src="assets/people/{esc(m["image"])}" alt="{esc(m["name"])}" loading="lazy"><span>{esc(m["name"])}</span></div>' for m in people["mascots"]),
        "{{NEWS}}": news_items(news),
        "{{NEWS_RECENT}}": news_items(news[:3]),
        "{{RESOURCES}}": resource_groups(resources),
        "{{FUNDERS}}": "\n".join(f'<img src="assets/funders/{f}" alt="{esc(n)}" loading="lazy">' for f, n in FUNDERS),
        "{{YEAR}}": str(YEAR),
    }
    for frag in sorted((SRC / "pages").glob("*.html")):
        text = frag.read_text(encoding="utf-8")
        m = re.match(r"<!--\s*title:\s*(.*?)\s*\|\s*desc:\s*(.*?)\s*-->", text)
        title, desc = (m.group(1), m.group(2)) if m else (frag.stem.title(), "")
        body = text[m.end():] if m else text
        for k, v in fills.items(): body = body.replace(k, v)
        (OUT / frag.name).write_text(head(title, desc, frag.name) + body + footer(), encoding="utf-8")
        print("wrote", frag.name)
    (OUT / "publications.html").write_text(
        head("Publications", "Publications from the Gawad Lab at Stanford Medicine: single-cell genomics, cell-free DNA, and pediatric leukemia.", "publications.html")
        + publications_page(pubs) + footer(), encoding="utf-8")
    print("wrote publications.html")
    pages = [p.name for p in OUT.glob("*.html")]
    (OUT / "sitemap.xml").write_text('<?xml version="1.0" encoding="UTF-8"?>\n<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
        + "".join(f"  <url><loc>{SITE_URL}/{'' if p=='index.html' else p}</loc></url>\n" for p in sorted(pages)) + "</urlset>\n")
    (OUT / "robots.txt").write_text(f"User-agent: *\nAllow: /\nSitemap: {SITE_URL}/sitemap.xml\n")
    (OUT / "CNAME").write_text("gawadlab.org\n")
    (OUT / ".nojekyll").write_text("")
    print(f"done: {len(pages)} pages, {len(pubs)} publications")

if __name__ == "__main__":
    build()
