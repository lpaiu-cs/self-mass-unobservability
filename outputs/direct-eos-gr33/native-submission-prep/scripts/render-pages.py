"""Render selected PDF pages to PNG for visual QA. Usage: python render-pages.py <pdf> <first> <last> [scale]"""
import sys
from pathlib import Path
import pypdfium2 as pdfium
pdf = pdfium.PdfDocument(sys.argv[1]); first, last = int(sys.argv[2]), int(sys.argv[3]); scale = float(sys.argv[4]) if len(sys.argv) > 4 else 1.5
out = Path(__file__).resolve().parent/'pages'; out.mkdir(exist_ok=True)
print('pages', len(pdf))
for i in range(first, last + 1):
    img = pdf[i - 1].render(scale=scale).to_pil(); p = out/f'page-{i:02d}.png'; img.save(p); print(p, img.size)
