"""Print the Ransom et al. 2014 (arXiv:1401.0535) lines that give e sin/cos omega, omega and ascending-node times."""
import re, sys
import pypdfium2 as pdfium
d = pdfium.PdfDocument(sys.argv[1])
pat = re.compile(r'sin|cos|97\.6|95\.6|scending|eriastron|ω')
for i in range(len(d)):
    t = d[i].get_textpage().get_text_range()
    if 'sin' in t and ('97.6' in t or '95.6' in t):
        for line in t.splitlines():
            if pat.search(line): print(i + 1, '|', line.strip()[:170])
