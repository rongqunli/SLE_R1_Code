import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
from pathlib import Path
import hashlib, json, subprocess
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader

root = Path(__file__).resolve().parents[1]
out = root / "work"
out.mkdir(parents=True, exist_ok=True)
records = root / "figures"
source = Path(__file__)
W, H = (130 / 25.4 * 72, 125 / 25.4 * 72)
margin, gap = (12.0, 8.0)
col = (W - 2 * margin - 2 * gap) / 3
xs = [margin + i * (col + gap) for i in range(3)]
cx = [x + col / 2 for x in xs]
pdfmetrics.registerFont(TTFont("Arial", "C:/Windows/Fonts/arial.ttf"))
pdfmetrics.registerFont(TTFont("Arial-Bold", "C:/Windows/Fonts/arialbd.ttf"))
INK = "#273746"
GREY = "#687784"
BORDER = "#B6C0C7"
PALE = "#F7F9FA"
MAGENTA = "#B53C84"
GREEN = "#729436"
CYAN = "#268FA2"
c = canvas.Canvas(str(out / "Figure0.pdf"), pagesize=(W, H), pageCompression=1)
c.setTitle("Study overview of the SLE immune-program analyses")
c.setAuthor("SLE R1 analysis team")
c.setSubject("Schematic workflow corresponding to the accepted R1 manuscript")
c.setFillColor(white)
c.rect(0, 0, W, H, stroke=0, fill=1)
elements = []


def line(x1, y1, x2, y2, color=GREY, width=0.8):
    c.setStrokeColor(HexColor(color))
    c.setLineWidth(width)
    c.line(x1, H - y1, x2, H - y2)


def arrow(x, y1, y2, color=GREY):
    line(x, y1, x, y2 - 3.8, color)
    p = c.beginPath()
    p.moveTo(x, H - y2)
    p.lineTo(x - 2.8, H - y2 + 4.2)
    p.lineTo(x + 2.8, H - y2 + 4.2)
    p.close()
    c.setFillColor(HexColor(color))
    c.drawPath(p, stroke=0, fill=1)


def text(x, y, string, size=8.5, bold=False, color=INK, maxwidth=None):
    font = "Arial-Bold" if bold else "Arial"
    width = pdfmetrics.stringWidth(string, font, size)
    if maxwidth is not None:
        assert width <= maxwidth, (string, width, maxwidth)
    c.setFont(font, size)
    c.setFillColor(HexColor(color))
    c.drawCentredString(x, H - y, string)
    elements.append(
        {
            "text": string,
            "font_size_pt": size,
            "center_x": x,
            "baseline_from_top": y,
            "width_pt": width,
        }
    )


def box(x, y, w, h, stroke=BORDER, fill=PALE):
    c.setStrokeColor(HexColor(stroke))
    c.setFillColor(HexColor(fill))
    c.setLineWidth(0.8)
    c.roundRect(x, H - y - h, w, h, 4.5, stroke=1, fill=1)


box(55, 10, W - 110, 48)
text(W / 2, 24, "Pediatric longitudinal discovery", 10, bold=True)
text(W / 2, 37, "GSE65391 | 131 patients, 897 visits", 8.4)
text(W / 2, 49, "Within-patient coexpression network", 8.4)
line(W / 2, 58, W / 2, 68)
line(cx[0], 68, cx[2], 68)
for x in cx:
    arrow(x, 68, 80)
for x, color, fill in zip(xs, [MAGENTA, GREEN, CYAN], ["#FCF5F9", "#F7FAF1", "#F2FAFB"]):
    box(x, 80, col, 43, color, fill)
text(cx[0], 96, "IFN / antiviral", 9.1, True, MAGENTA, col - 10)
text(cx[0], 111, "Magenta module", 8.1, maxwidth=col - 10)
text(cx[1], 96, "B-cell / BCR", 9.1, True, GREEN, col - 10)
text(cx[1], 111, "Greenyellow module", 8.1, maxwidth=col - 10)
text(cx[2], 93, "Immature-neutrophil-", 8.2, True, CYAN, col - 8)
text(cx[2], 103, "associated program", 8.2, True, CYAN, col - 8)
text(cx[2], 115, "Cyan module", 8.1, maxwidth=col - 10)
for x in cx:
    line(x, 123, x, 133)
line(cx[0], 133, cx[2], 133)
adult_gap = 12.0
adult_w = (W - 2 * margin - adult_gap) / 2
adult_x = [margin, margin + adult_w + adult_gap]
adult_c = [x + adult_w / 2 for x in adult_x]
for x in adult_c:
    arrow(x, 133, 152)
    box(x - adult_w / 2, 152, adult_w, 60)
text(adult_c[0], 168, "Adult network comparison", 9.1, True, maxwidth=adult_w - 12)
text(adult_c[0], 181, "GSE88884 | 1,319 discovery participants", 7.8, maxwidth=adult_w - 12)
text(adult_c[0], 194, "Gene overlap, network preservation", 8.1, maxwidth=adult_w - 12)
text(adult_c[0], 205, "and module-membership concordance", 8.1, maxwidth=adult_w - 12)
text(adult_c[1], 168, "Continuous clinical coordinates", 8.7, True, maxwidth=adult_w - 12)
text(adult_c[1], 181, "Joint IFN, B-cell and neutrophil scores", 8.0, maxwidth=adult_w - 12)
text(adult_c[1], 194, "Clinical associations and held-out evaluation", 7.5, maxwidth=adult_w - 12)
text(adult_c[1], 205, "Reference-group and IFN6 benchmarks", 8.0, maxwidth=adult_w - 12)
for x in adult_c:
    line(x, 212, x, 225)
line(cx[0], 225, cx[2], 225)
for x in cx:
    arrow(x, 225, 241)
for x, stroke, fill in zip(xs, [MAGENTA, GREEN, BORDER], ["#FCF5F9", "#F7FAF1", PALE]):
    box(x, 241, col, 98, stroke, fill)
text(cx[0], 257, "Cell-resolved IFN", 9.0, True, MAGENTA, col - 10)
text(cx[0], 272, "GSE135779", 8.2, maxwidth=col - 10)
text(cx[0], 284, "44 pediatric donors", 8.2, maxwidth=col - 10)
line(xs[0] + 10, 291, xs[0] + col - 10, 291, BORDER, 0.5)
text(cx[0], 304, "Cell-type effects", 8.1, maxwidth=col - 10)
text(cx[0], 316, "Donor-mode factorization", 7.9, maxwidth=col - 10)
text(cx[0], 328, "TF-expression associations", 7.9, maxwidth=col - 10)
text(cx[1], 254, "Trial molecular", 9.0, True, GREEN, col - 10)
text(cx[1], 265, "subsets", 9.0, True, GREEN, col - 10)
text(cx[1], 279, "ILLUMINATE-1 / -2", 8.0, maxwidth=col - 10)
text(cx[1], 291, "309 / 147 participants", 8.2, maxwidth=col - 10)
line(xs[1] + 10, 298, xs[1] + col - 10, 298, BORDER, 0.5)
text(cx[1], 311, "Tabalumab ANCOVA", 8.0, maxwidth=col - 10)
text(cx[1], 323, "B-cell abundance", 8.0, maxwidth=col - 10)
text(cx[1], 333, "sensitivity analysis", 8.0, maxwidth=col - 10)
text(cx[2], 254, "Adult PBMC", 9.0, True, maxwidth=col - 10)
text(cx[2], 265, "associations", 9.0, True, maxwidth=col - 10)
text(cx[2], 279, "GSE121239", 8.2, maxwidth=col - 10)
text(cx[2], 291, "65 patients, 292 visits", 8.2, maxwidth=col - 10)
line(xs[2] + 10, 298, xs[2] + col - 10, 298, BORDER, 0.5)
text(cx[2], 311, "GZMK-state / cytotoxicity", 7.8, maxwidth=col - 10)
text(cx[2], 323, "IFN-adjusted", 8.0, maxwidth=col - 10)
text(cx[2], 333, "associations", 8.0, maxwidth=col - 10)
c.showPage()
c.save()
subprocess.run(
    [
        "pdftoppm",
        "-png",
        "-r",
        "600",
        "-singlefile",
        str(out / "Figure0.pdf"),
        str(out / "Figure0"),
    ],
    check=True,
    capture_output=True,
)
subprocess.run(
    [
        "pdftoppm",
        "-png",
        "-r",
        "200",
        "-singlefile",
        str(out / "Figure0.pdf"),
        str(records / "Figure0_preview"),
    ],
    check=True,
    capture_output=True,
)
extract = "\n".join((page.extract_text() for page in PdfReader(out / "Figure0.pdf").pages))
for forbidden in [
    "trajectory",
    "trajectories",
    "milo",
    "regulon",
    "knockout",
    "perturbation",
    "b-axis",
    "share",
]:
    assert forbidden not in extract.lower(), forbidden
assert "897" in extract and "131" in extract and ("309 / 147" in extract) and ("292" in extract)
(records / "Figure0_extracted_text.txt").write_text(extract, encoding="utf-8")
manifest = {
    "source": str(source),
    "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "kind": "schematic workflow; no data plots",
    "size_mm": [130, 125],
    "png_dpi": 600,
    "nodes": elements,
    "files": {
        name: hashlib.sha256((out / name).read_bytes()).hexdigest()
        for name in ["Figure0.pdf", "Figure0.png"]
    },
}
(records / "Figure0_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(
    json.dumps(
        {
            "outputs": [str(out / "Figure0.pdf"), str(out / "Figure0.png")],
            "size_mm": [130, 125],
            "minimum_font_size_pt": min((e["font_size_pt"] for e in elements)),
        },
        indent=2,
    )
)
