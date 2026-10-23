"""Render the authored corrected result report without network dependencies."""
from pathlib import Path
import html
import re
import shutil

ROOT = Path(__file__).resolve().parent.parent


def inline(text):
    text = html.escape(text)
    text = re.sub(r'`([^`]+)`', r'<code>\1</code>', text)
    return re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', text)


def body(markdown):
    output, paragraph, table, bullets = [], [], [], []
    def flush():
        if paragraph:
            output.append('<p>' + inline(' '.join(paragraph)) + '</p>')
            paragraph.clear()
        if table:
            header, *rows = table
            output.append('<div class="table-wrap"><table><thead><tr>' + ''.join('<th>'+inline(c)+'</th>' for c in header) + '</tr></thead><tbody>')
            for row in rows:
                output.append('<tr>' + ''.join('<td>'+inline(c)+'</td>' for c in row) + '</tr>')
            output.append('</tbody></table></div>')
            table.clear()
        if bullets:
            output.append('<ul>' + ''.join('<li>'+inline(c)+'</li>' for c in bullets) + '</ul>')
            bullets.clear()
    for line in markdown.splitlines():
        if line.startswith('|'):
            if paragraph or bullets:
                flush()
            cells = [c.strip() for c in line.strip().strip('|').split('|')]
            if not all(re.match(r'^:?-+:?$', c) for c in cells):
                table.append(cells)
        elif line.startswith('- '):
            if paragraph or table:
                flush()
            bullets.append(line[2:])
        elif line.startswith('#'):
            flush()
            prefix, _, content = line.partition(' ')
            level = min(len(prefix), 6)
            output.append(f'<h{level}>'+inline(content)+f'</h{level}>')
        elif not line.strip():
            flush()
        else:
            if table or bullets:
                flush()
            paragraph.append(line)
    flush()
    return '\n'.join(output)


def build():
    target = ROOT / 'site/communication'
    # Historical report is kept in Git and immutable v3.5 tags. Its outdated
    # biological selection is deliberately not another active report here.
    source = ROOT / 'docs/DEG_ONLY_RESULTS.md'
    shutil.copyfile(source, target / 'RESULTS_AND_INTERPRETATION.md')
    document = '''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Corrected DEG-only CA1–inhibitory results</title><link rel="stylesheet" href="styles.css"><style>.report{max-width:1200px;margin:2rem auto;padding:2.2rem;background:white;border:1px solid #d7e2e8;border-radius:12px}.report h2{margin-top:2.2rem}.report h3{margin-top:1.5rem}.report table{font-size:.85rem}.report nav{display:flex;gap:1rem;flex-wrap:wrap;margin-bottom:1.8rem}.report code{font-size:.87em;overflow-wrap:anywhere;background:#edf3f7;padding:.1rem .25rem}.report ul{line-height:1.8;padding-left:1.3rem}.report .table-wrap{max-height:none}@media(max-width:700px){.report{margin:.5rem;padding:1rem}}</style></head><body><main class="report"><nav><a href="index.html">DEG-only interactive atlas</a><a href="data/deg_only/index.html">Corrected downloads and exact DEGs</a><a href="RESULTS_AND_INTERPRETATION.md" download>Download report</a></nav><div class="notice">Corrected selection · atlas v3.6.0 · 23 October 2026. Both matching endpoint genes must be DEGs. Native LIANA scores retain their original full-background meaning; no scientific fit was rerun.</div>'''
    (target / 'RESULTS.html').write_text(document + body(source.read_text()) + '</main></body></html>')


if __name__ == '__main__':
    build()
