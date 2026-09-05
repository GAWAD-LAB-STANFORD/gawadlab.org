# Gawad Lab website (gawadlab.org)

A modern, dependency-free static rebuild of gawadlab.org. Plain HTML, CSS, and a
few dozen lines of JavaScript. No WordPress, no plugins, nothing to patch.

## Layout

```
build.py                    builds the site into site/
scripts/update_publications.py   pulls new papers from PubMed and rebuilds
data/publications.json      every paper (generated from PubMed, hand-editable)
data/people.json            current members, alumni, mascots
data/news.json              news items (newest first)
data/resources.json         links on the Resources page
src/pages/*.html            page templates (index, research, people, join, news, resources)
src/css/style.css           the stylesheet
src/js/main.js              mobile menu, publication search, author toggles
src/assets/                 logo, hero, headshots (assets/people), funder logos (assets/funders)
site/                       the built site: upload this folder as-is
```

## Everyday edits

* **New paper**: `python3 scripts/update_publications.py`. It adds anything new
  from PubMed and rebuilds. Open `data/publications.json` to set
  `"highlight": true` on papers you want in the "Selected" filter, and delete
  a preprint once the journal version appears.
* **New lab member**: add an entry to `data/people.json` and drop a square
  headshot (about 400x400 px) into `src/assets/people/`. Move departing members
  to `alumni` with their years and where they went.
* **News**: add an item at the top of `data/news.json`.
* **Research text**: edit `src/pages/research.html` (and the five cards in
  `src/pages/index.html`).
* Then run `python3 build.py` and upload `site/`.

## Publishing options

1. **GitHub Pages (recommended, free)**: push this folder to a repo under
   github.com/GAWAD-LAB-STANFORD, enable Pages from the `site/` folder (or a
   `docs/` folder), and point gawadlab.org's DNS at GitHub Pages. HTTPS is
   automatic.
2. **Keep GoDaddy hosting**: upload the contents of `site/` to the web root
   (replacing WordPress). Static files need no PHP or database.
3. **Netlify / Cloudflare Pages**: drag-and-drop `site/`, then add the custom
   domain.

## Before going live

* Verify the roster in `data/people.json`; it was carried over from the 2022
  site.
* Verify the funder logos on the home page reflect current support.
* Read the research program text on `research.html`; it was written from
  published work and public lab materials and should be checked for anything
  you would rather not state publicly.
* Preview locally with `python3 -m http.server 8080 --directory site` and open
  http://localhost:8080.
