"""Build the compact Pages entry from the editable poster package.

Run with 'body', commit the package, then run with that commit SHA to build the
Pages entry. The static asset pin prevents a stale editor from loading new data.
"""
import base64
import hashlib
import json
from html import escape
from html.parser import HTMLParser
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "analysis" / "neuroday26"
ENTRY = ROOT / "site" / "neuroday26" / "index.html"


class Resources(HTMLParser):
    def __init__(self):
        super().__init__()
        self.css = []
        self.scripts = []
        self.title = ""
        self.in_title = False

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "title":
            self.in_title = True
        elif tag == "link" and a.get("rel") == "stylesheet":
            self.css.append(a["href"])
        elif tag == "script":
            self.scripts.append(a["src"])

    def handle_endtag(self, tag):
        if tag == "title":
            self.in_title = False

    def handle_data(self, data):
        if self.in_title:
            self.title += data


def main():
    if len(sys.argv) != 2:
        sys.exit("Expected 'body' or the full asset commit SHA.")
    rev = sys.argv[1]
    source = (PACKAGE / "index.html").read_text()
    resources = Resources()
    resources.feed(source)
    body = re.search(r"<body>(.*?)</body>", source, re.DOTALL | re.IGNORECASE).group(1)
    body = re.sub(r"<script\b[^>]*>.*?</script>", "", body, flags=re.DOTALL | re.IGNORECASE)
    bootstrap = PACKAGE / "atlas-body.js"
    if rev == "body":
        # No HTML base tag: local SVG fragment references retain the page origin.
        bootstrap.write_text(
            'window.POSTER_ASSET_BASE = new URL(".", document.currentScript.src).href;\n'
            'const posterUI = document.createElement("template");\n'
            + "posterUI.innerHTML = " + json.dumps(body, ensure_ascii=False) + ";\n"
            # Individual figure files are served directly by GitHub. The CDN also
            # serves editor scripts and fonts but can decline an uncached figure
            # from a large repository even when the figure itself is small.
            'const packageURL = new URL(window.POSTER_ASSET_BASE);\n'
            'const parts = packageURL.pathname.split("/");\n'
            'const revision = parts[3]?.split("@");\n'
            'const isPinned = packageURL.hostname === "cdn.jsdelivr.net" && parts[1] === "gh" && revision?.length === 2 && /^[a-f0-9]{40}$/.test(revision[1]);\n'
            'const figureBase = isPinned ? new URL([parts[2], revision[0], revision[1], ...parts.slice(4)].join("/"), "https://raw.githubusercontent.com/").href : window.POSTER_ASSET_BASE;\n'
            'for (const link of posterUI.content.querySelectorAll("a[href]")) {\n'
            '  const href = link.getAttribute("href");\n'
            '  const base = link.classList.contains("tree-source") ? figureBase : window.POSTER_ASSET_BASE;\n'
            '  if (href && !href.startsWith("#")) link.href = new URL(href, base).href;\n'
            '}\ndocument.currentScript.before(posterUI.content);\n'
        )
        return
    if not re.fullmatch(r"[0-9a-f]{40}", rev):
        sys.exit("Expected a complete hexadecimal commit SHA.")
    base = f"https://cdn.jsdelivr.net/gh/Ali-Mujahid/satb2_snRNA@{rev}/analysis/neuroday26/"
    source_base = f"https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/{rev}/analysis/neuroday26/"

    def asset(name):
        path = (PACKAGE / name).resolve()
        relative = path.relative_to(PACKAGE).as_posix()
        integrity = base64.b64encode(hashlib.sha384(path.read_bytes()).digest()).decode()
        return f'"{base}{relative}" integrity="sha384-{integrity}" crossorigin="anonymous"'

    styles = "\n".join(f'<link rel="stylesheet" href={asset(path)}>' for path in resources.css)
    scripts = "\n".join(f'<script defer src={asset(path)}></script>' for path in resources.scripts)
    html = (
        '<!doctype html>\n<html lang="en"><head><meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width,initial-scale=1">\n'
        f'<title>{escape(resources.title)}</title>\n{styles}\n{scripts}\n</head>\n'
        f'<body><script src={asset(bootstrap.name)}></script>\n'
        '<noscript>Enable JavaScript to edit this poster. Download the original '
        f'<a href="{source_base}poster.svg">SVG</a> or <a href="{source_base}poster.pdf">PDF</a>.</noscript>\n'
        '</body></html>\n'
    )
    ENTRY.parent.mkdir(parents=True, exist_ok=True)
    ENTRY.write_text(html)


if __name__ == "__main__":
    main()
