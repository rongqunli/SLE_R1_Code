import sys

sys.dont_write_bytecode = True
import argparse
import json
import math
from pathlib import Path

import pypdfium2 as pdfium
from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject


CODE_ROOT = Path(__file__).resolve().parents[2]


def output_path(value, inputs=()):
    path = Path(value).expanduser().resolve()
    if path == CODE_ROOT or path.is_relative_to(CODE_ROOT):
        raise ValueError("Generated files must be outside the code repository")
    if path in {Path(item).expanduser().resolve() for item in inputs}:
        raise ValueError("An output must not overwrite an input PDF")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def rasterize(source, dpi, opaque=False):
    if not math.isfinite(dpi) or dpi <= 0:
        raise ValueError("DPI must be positive and finite")
    with pdfium.PdfDocument(str(source)) as document:
        if len(document) != 1:
            raise ValueError("Expected a one-page PDF: " + str(source))
        page = document[0]
        size = page.get_size()
        bitmap = page.render(
            scale=dpi / 72,
            fill_color=(255, 255, 255, 255) if opaque else (0, 0, 0, 0),
        )
        image = bitmap.to_pil().convert("RGB" if opaque else "RGBA").copy()
        bitmap.close()
        page.close()
    return image, size


def render(source, destination, dpi=600, opaque=False):
    target = output_path(destination, [source])
    image, _ = rasterize(source, dpi, opaque)
    image.save(target, dpi=(dpi, dpi))
    return target


def single_page(source):
    reader = PdfReader(source)
    if len(reader.pages) != 1:
        raise ValueError("Expected a one-page PDF: " + str(source))
    page = reader.pages[0]
    if page.rotation:
        page.transfer_rotation_to_content()
    return page


def save_pdf(page, destination, inputs):
    target = output_path(destination, inputs)
    writer = PdfWriter()
    writer.add_page(page)
    writer.pages[0].compress_content_streams()
    with target.open("wb") as handle:
        writer.write(handle)
    return target


def crop_pair(hd, sle, destination, dpi=600, padding=3):
    if not math.isfinite(padding) or padding < 0:
        raise ValueError("Padding must be nonnegative and finite")
    sources = {"HD": Path(hd), "SLE": Path(sle)}
    records = {}
    for group, source in sources.items():
        image, size = rasterize(source, dpi)
        box = image.getchannel("A").getbbox()
        if box is None:
            raise ValueError("The PDF has no visible content: " + str(source))
        records[group] = {"page_size_pt": size, "alpha_bbox_pixels": box}
    width, height = records["HD"]["page_size_pt"]
    if any(abs(a - b) > 1e-6 for a, b in zip(records["SLE"]["page_size_pt"], (width, height))):
        raise ValueError("The HD and SLE uncropped pages must have identical dimensions")
    boxes = [record["alpha_bbox_pixels"] for record in records.values()]
    scale = dpi / 72
    left = max(0, min(box[0] for box in boxes) / scale - padding)
    top = max(0, min(box[1] for box in boxes) / scale - padding)
    right = min(width, max(box[2] for box in boxes) / scale + padding)
    bottom = min(height, max(box[3] for box in boxes) / scale + padding)
    output_size = (right - left, bottom - top)
    for group, source in sources.items():
        page = single_page(source)
        if any(abs(float(value)) > 1e-6 for value in page.mediabox.lower_left):
            raise ValueError("Expected an uncropped PDF page with zero origin")
        page.add_transformation(Transformation().translate(tx=-left, ty=-(height - bottom)))
        page.mediabox = RectangleObject((0, 0, *output_size))
        page.cropbox = RectangleObject((0, 0, *output_size))
        for key in ("/TrimBox", "/BleedBox", "/ArtBox"):
            if key in page:
                del page[key]
        target = Path(destination) / f"Figure5B_Circular_{group}_for_AI.pdf"
        target = save_pdf(page, target, sources.values())
        png = render(target, target.with_suffix(".png"), dpi)
        with Image.open(png) as image:
            alpha = image.getchannel("A")
            box = alpha.getbbox()
            if box is None:
                raise ValueError("The cropped PDF has no visible content")
            if any(alpha.getpixel(point) != 0 for point in ((0, 0), (image.width - 1, 0), (0, image.height - 1), (image.width - 1, image.height - 1))):
                raise ValueError("Cropped single-group output lost corner transparency")
            records[group].update({"pdf": str(target), "png": str(png), "png_pixels": image.size, "content_bbox_pixels": box})
    record = {
        "dpi": dpi,
        "padding_pt": padding,
        "common_crop_top_origin_pt": [left, top, right, bottom],
        "common_page_size_pt": output_size,
        "groups": records,
    }
    manifest = output_path(Path(destination) / "Figure5B_crop_geometry.json", sources.values())
    manifest.write_text(json.dumps(record, indent=2), encoding="utf-8")
    return record


def assemble(base, panel_a, panel_b, destination, top_height=118.8, overlay=None):
    page = single_page(base)
    width, height = float(page.mediabox.width), float(page.mediabox.height)
    if not math.isfinite(top_height) or not 0 < top_height < height:
        raise ValueError("The top-row height must fit inside the base page")
    for index, source in enumerate((panel_a, panel_b)):
        panel = single_page(source)
        box = panel.cropbox
        panel_width, panel_height = float(box.width), float(box.height)
        factor = min(width / 2 / panel_width, top_height / panel_height)
        x = index * width / 2 + (width / 2 - panel_width * factor) / 2
        y = height - top_height + (top_height - panel_height * factor) / 2
        matrix = Transformation().translate(tx=-float(box.left), ty=-float(box.bottom)).scale(factor).translate(tx=x, ty=y)
        page.merge_transformed_page(panel, matrix, over=True, expand=False)
    inputs = [base, panel_a, panel_b]
    if overlay:
        annotation = single_page(overlay)
        if any(abs(float(a) - float(b)) > 1e-6 for a, b in zip(annotation.mediabox, page.mediabox)):
            raise ValueError("The annotation overlay must have the same page box as the base")
        page.merge_page(annotation, over=True, expand=False)
        inputs.append(overlay)
    return save_pdf(page, destination, inputs)


def main():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    renderer = commands.add_parser("render")
    renderer.add_argument("input")
    renderer.add_argument("output")
    renderer.add_argument("--dpi", type=float, default=600)
    renderer.add_argument("--opaque", action="store_true")
    cropper = commands.add_parser("crop-pair")
    cropper.add_argument("hd")
    cropper.add_argument("sle")
    cropper.add_argument("output_directory")
    cropper.add_argument("--dpi", type=float, default=600)
    cropper.add_argument("--padding", type=float, default=3)
    assembler = commands.add_parser("assemble")
    assembler.add_argument("base")
    assembler.add_argument("panel_a")
    assembler.add_argument("panel_b")
    assembler.add_argument("output")
    assembler.add_argument("--top-height", type=float, default=118.8)
    assembler.add_argument("--overlay")
    args = parser.parse_args()
    if args.command == "render":
        render(args.input, args.output, args.dpi, args.opaque)
    elif args.command == "crop-pair":
        crop_pair(args.hd, args.sle, args.output_directory, args.dpi, args.padding)
    else:
        assemble(args.base, args.panel_a, args.panel_b, args.output, args.top_height, args.overlay)


if __name__ == "__main__":
    main()
