#!/usr/bin/env python3
import argparse, datetime, sys
from pathlib import Path

import numpy as np
from tifffile import TiffFile, imwrite
from PIL import Image, ImageOps, ImageEnhance

# Optional: OpenCV CLAHE for better local contrast if available
try:
    import cv2
    HAVE_CV2 = True
except Exception:
    HAVE_CV2 = False

def get_tag(page, code):
    try:
        return page.tags[code].value
    except KeyError:
        return None

def pil_autocontrast_rgb(img, cutoff=0.5, gain=1.0):
    # img: HxWx3 uint8
    out = np.empty_like(img)
    for c in range(3):
        ch = Image.fromarray(img[..., c], mode="L")
        ch = ImageOps.autocontrast(ch, cutoff=cutoff)
        ch = ImageEnhance.Contrast(ch).enhance(gain)
        out[..., c] = np.array(ch, dtype=np.uint8)
    return out

def cv2_clahe_rgb(img, clip_limit=2.0, tile_grid_size=8):
    if not HAVE_CV2:
        raise RuntimeError("OpenCV not available; use --method autocontrast or install opencv-python")
    clahe = cv2.createCLAHE(clipLimit=float(clip_limit), tileGridSize=(tile_grid_size, tile_grid_size))
    out = np.empty_like(img)
    for c in range(3):
        out[..., c] = clahe.apply(img[..., c])
    return out

def main():
    ap = argparse.ArgumentParser(description="Enhance CytAssist/high-res TIFF contrast without losing key metadata.")
    ap.add_argument("input", help="Input TIFF path (CytAssist or high-res).")
    ap.add_argument("--output", help="Output TIFF path (default: <input>.contrast.tif)")
    ap.add_argument("--method", choices=["autocontrast", "clahe"], default="autocontrast",
                    help="Contrast method: autocontrast (Pillow) or clahe (OpenCV).")
    ap.add_argument("--autocontrast-cutoff", type=float, default=0.5, help="Percent cutoff per tail (autocontrast).")
    ap.add_argument("--autocontrast-gain", type=float, default=1.1, help="Global contrast gain after autocontrast.")
    ap.add_argument("--clahe-clip", type=float, default=2.0, help="CLAHE clip limit (higher = more contrast).")
    ap.add_argument("--clahe-tile", type=int, default=8, help="CLAHE tile grid size (N x N).")
    args = ap.parse_args()

    ipath = Path(args.input)
    if not ipath.exists():
        sys.exit(f"Input not found: {ipath}")

    opath = Path(args.output) if args.output else ipath.with_suffix(".contrast.tif")

    # Read first page; CytAssist and high-res are usually single-page RGB, uint8
    with TiffFile(str(ipath)) as tf:
        page = tf.pages[0]
        arr = page.asarray()
        # Basic sanity
        assert arr.ndim == 3 and arr.shape[2] in (3, 4), "Expect HxWx3 or HxWx4 image"
        if arr.shape[2] == 4:
            # Drop alpha if present (Space Ranger expects RGB)
            arr = arr[..., :3]
        if arr.dtype != np.uint8:
            # Scale/clip to uint8 conservatively
            arr = np.clip(arr, 0, 255).astype(np.uint8)

        # Capture key tags to preserve
        # TIFF tag IDs: 282 XResolution, 283 YResolution, 296 ResolutionUnit, 274 Orientation,
        # 305 Software, 306 DateTime, 270 ImageDescription
        xres = get_tag(page, 282)
        yres = get_tag(page, 283)
        resunit = get_tag(page, 296)
        orient = get_tag(page, 274)
        software = get_tag(page, 305)
        dt = get_tag(page, 306)
        desc = get_tag(page, 270)

    # Enhance
    if args.method == "autocontrast":
        enhanced = pil_autocontrast_rgb(arr, cutoff=args.autocontrast_cutoff, gain=args.autocontrast_gain)
    else:
        enhanced = cv2_clahe_rgb(arr, clip_limit=args.clahe_clip, tile_grid_size=args.clahe_tile)

    # Prepare extratags to carry through key metadata if present
    # extratag tuple: (code, dtype, count, value, writeonce)
    # dtype: 5=RATIONAL, 3=SHORT, 2=ASCII
    extratags = []
    def add_tag(code, dtype, value):
        if value is None:
            return
        # Normalize common forms
        if code in (282, 283):  # RATIONAL (X/Y Resolution)
            # tifffile accepts tuples like (num, den) or float
            if isinstance(value, (tuple, list)) and len(value) == 2:
                pass
            elif isinstance(value, (int, float)):
                value = (int(value*10000), 10000)
            else:
                # best-effort default 72 dpi if value weird
                value = (72, 1)
            extratags.append((code, 5, 1, value, False))
        elif code in (296, 274):  # SHORT
            try:
                val = int(value)
            except Exception:
                return
            extratags.append((code, 3, 1, val, False))
        elif code in (305, 306, 270):  # ASCII
            s = str(value)
            extratags.append((code, 2, len(s) + 1, s, False))

    add_tag(282, 5, xres)
    add_tag(283, 5, yres)
    add_tag(296, 3, resunit)
    add_tag(274, 3, orient)
    add_tag(305, 2, software if software else "contrast-enhanced by script")
    add_tag(306, 2, dt if dt else datetime.datetime.now().strftime("%Y:%m:%d %H:%M:%S"))
    add_tag(270, 2, desc)

    # Write: LZW, RGB, keep tags; avoid ImageJ metadata
    imwrite(
        str(opath),
        enhanced,
        photometric="rgb",
        compression="LZW",
        metadata=None,
        extratags=extratags
    )

    # Minimal report
    print(f"Wrote: {opath}")
    print(f"Size: {enhanced.shape[1]}x{enhanced.shape[0]}, dtype: {enhanced.dtype}, method: {args.method}")

if __name__ == "__main__":
    main()
