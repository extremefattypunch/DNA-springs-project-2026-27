#!/usr/bin/env python3
"""Crop a region of a gel image, upscale it and draw a labelled pixel grid.

Coordinates on the grid are ORIGINAL image pixels, so you can read band and
lane positions straight off the picture.

usage: grid.py IMAGE OUT.png [--box x y w h] [--scale 4] [--step 10] [--invert]
"""
import argparse
from PIL import Image, ImageDraw, ImageFont, ImageOps

ap = argparse.ArgumentParser()
ap.add_argument("image")
ap.add_argument("out")
ap.add_argument("--box", nargs=4, type=int, metavar=("X", "Y", "W", "H"))
ap.add_argument("--scale", type=float, default=4)
ap.add_argument("--step", type=int, default=10, help="grid spacing in original px")
ap.add_argument("--invert", action="store_true")
a = ap.parse_args()

im = Image.open(a.image).convert("RGB")
if a.box:
    x, y, w, h = a.box
else:
    x, y, w, h = 0, 0, im.width, im.height
crop = im.crop((x, y, x + w, y + h))
if a.invert:
    crop = ImageOps.invert(crop)
s = a.scale
big = crop.resize((int(w * s), int(h * s)), Image.LANCZOS)
pad_l, pad_t = 44, 22
canvas = Image.new("RGB", (big.width + pad_l, big.height + pad_t), (255, 255, 255))
canvas.paste(big, (pad_l, pad_t))
d = ImageDraw.Draw(canvas)
try:
    font = ImageFont.truetype("/usr/share/fonts/TTF/DejaVuSans.ttf", 11)
except OSError:
    font = ImageFont.load_default()

start_x = (x // a.step + 1) * a.step if x % a.step else x
for gx in range(start_x, x + w + 1, a.step):
    px = pad_l + (gx - x) * s
    major = gx % (a.step * 5) == 0
    d.line([(px, pad_t), (px, canvas.height)], fill=(255, 60, 60) if major else (255, 170, 0), width=1)
    if major or s * a.step >= 40:
        d.text((px + 2, 4), str(gx), fill=(200, 0, 0), font=font)
start_y = (y // a.step + 1) * a.step if y % a.step else y
for gy in range(start_y, y + h + 1, a.step):
    py = pad_t + (gy - y) * s
    major = gy % (a.step * 5) == 0
    d.line([(pad_l, py), (canvas.width, py)], fill=(255, 60, 60) if major else (0, 190, 255), width=1)
    if major or s * a.step >= 40:
        d.text((2, py - 6), str(gy), fill=(0, 0, 200), font=font)
canvas.save(a.out)
print(f"wrote {a.out} ({canvas.width}x{canvas.height}); grid every {a.step}px of the original")
