"""Give changed frontend files distinct URLs while retaining the pinned data tag."""
from pathlib import Path
import hashlib,re

ROOT=Path(__file__).resolve().parent.parent
SITE=ROOT/'site'
for page in SITE.glob('*.html'):
    def replace(match):
        attr,path=match.group(1,2)
        source=page.parent/path
        if not source.is_file():return match.group(0)
        version=hashlib.sha256(source.read_bytes()).hexdigest()[:12]
        return f'{attr}="{path}?v={version}"'
    page.write_text(re.sub(r'(src|href)="([^"?#]+\.(?:js|css))(?:\?v=[0-9a-f]+)?"',replace,page.read_text()))
