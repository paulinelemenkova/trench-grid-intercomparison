#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compress_figures.py
===================
Make every PNG in a folder smaller than a size cap (default 2 MB), without
destroying quality, by downsampling only the oversized ones to a sensible print
resolution and re-saving with maximum PNG compression. Files already under the
cap are copied through unchanged.

NON-DESTRUCTIVE: originals are never touched. Results go to a sibling folder
(default: <figures>/compressed/). Check them, then swap them in.

Usage
-----
  python3 compress_figures.py                         # uses the path below
  python3 compress_figures.py /path/to/figures        # any folder
  python3 compress_figures.py /path/to/figures 1.5    # cap at 1.5 MB

Why downsampling (not JPEG): these are colour relief/map panels with text
labels. Downsampling to ~300 dpi at journal column width keeps them crisp in
print and slashes file size, while staying PNG so your \\includegraphics names
do not change. Line-art figures that are already small pass through untouched.
"""

import sys, os, shutil
from PIL import Image

# ---- settings ----
DEFAULT_DIR = "figures"   # relative default; pass your own path as argument 1
CAP_MB      = 1.5        # size cap per file (MB)
START_MAXPX = 2600       # first try: cap the longer side to this many pixels
FLOOR_MAXPX = 1800       # never shrink the longer side below this
STEP        = 0.85       # shrink factor per iteration if still too big
# A 190 mm (double-column) figure at 300 dpi is ~2244 px wide,
# so START_MAXPX=2600 preserves >=300 dpi; the 1800 floor is ~240 dpi at that
# width, reached only if a file cannot otherwise fit the cap.


def target_bytes(cap_mb):
    return int(cap_mb * 1024 * 1024)


def save_png(img, path):
    # optimize + max zlib compression; strip metadata to shave a little more
    img.save(path, format="PNG", optimize=True, compress_level=9)


def compress_one(src, dst, cap):
    size = os.path.getsize(src)
    if size <= cap:
        shutil.copy2(src, dst)
        return size, size, "copied (already under cap)"

    im = Image.open(src)
    im.load()
    if im.mode == "P":
        im = im.convert("RGBA")
    if im.mode == "RGBA":
        # flatten onto white so quantization/PNG stays clean
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[-1])
        im = bg
    elif im.mode != "RGB":
        im = im.convert("RGB")

    # 1) recompress at full resolution, no quality loss
    save_png(im, dst)
    if os.path.getsize(dst) <= cap:
        return size, os.path.getsize(dst), "recompressed, no resize"

    # 2) ordered ladder, least-degrading first. Each PASS steps the longer side
    #    from START down to FLOOR and returns the HIGHEST resolution that fits the
    #    cap; passes escalate only when a whole pass fails. So a dense relief map
    #    keeps ~300 dpi in 256 colours rather than dropping to a low-resolution or
    #    64-colour version. Colour depth is reduced only when full colour cannot
    #    fit at any resolution down to the floor.
    best = [os.path.getsize(dst), "recompressed full-res"]

    def resolutions():
        px = START_MAXPX
        while True:
            yield px
            if px <= FLOOR_MAXPX:
                return
            px = int(px * STEP)

    def write_variant(resized, colors):
        if colors is None:
            save_png(resized, dst)
        else:
            resized.quantize(colors=colors, method=Image.FASTOCTREE).save(
                dst, format="PNG", optimize=True, compress_level=9)
        return os.path.getsize(dst)

    for colors, label in ((None, "full colour"), (256, "256-colour palette"),
                          (128, "128-colour palette"), (64, "64-colour palette")):
        for px in resolutions():
            w, h = im.size
            scale = min(1.0, px / float(max(w, h)))
            rw, rh = max(1, int(round(w * scale))), max(1, int(round(h * scale)))
            resized = im.resize((rw, rh), Image.LANCZOS) if scale < 1.0 else im
            out = write_variant(resized, colors)
            if out <= cap:
                return size, out, f"resized to {rw}x{rh}, {label}"
            if out < best[0]:
                best = [out, f"{rw}x{rh}, {label}"]
    # nothing fit anywhere: dst holds the last (smallest) variant tried
    return size, best[0], best[1] + "  [STILL OVER CAP]"


def main():
    src_dir = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_DIR
    cap = target_bytes(float(sys.argv[2]) if len(sys.argv) > 2 else CAP_MB)
    if not os.path.isdir(src_dir):
        sys.exit(f"not a folder: {src_dir}")

    out_dir = os.path.join(src_dir, "compressed")
    os.makedirs(out_dir, exist_ok=True)

    pngs = sorted(f for f in os.listdir(src_dir) if f.lower().endswith(".png"))
    if not pngs:
        sys.exit(f"no PNG files in {src_dir}")

    print(f"{'file':<34}{'before':>10}{'after':>10}   note")
    print("-" * 78)
    over = []
    for f in pngs:
        b, a, note = compress_one(os.path.join(src_dir, f),
                                   os.path.join(out_dir, f), cap)
        print(f"{f:<34}{b/1e6:>8.2f}MB{a/1e6:>8.2f}MB   {note}")
        if a > cap:
            over.append(f)

    print("-" * 78)
    print(f"wrote {len(pngs)} files to {out_dir}")
    if over:
        print("STILL over cap (need a stronger step, see notes below):", over)
    else:
        print("all files are now under the cap.")


if __name__ == "__main__":
    main()
