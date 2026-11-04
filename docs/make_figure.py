# Builds docs/compare_example.png for the README: human alpha-globin,
# beta-globin and myoglobin compared with the "compare" command.
# Run from this folder:  ChimeraX --exit make_figure.py
import os
from chimerax.core.commands import run
from PIL import Image, ImageDraw, ImageFont

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)

for acc in ("P69905", "P68871", "P02144"):
    run(session, "alphafold fetch %s" % acc)
run(session, "rename #1 alpha-globin; rename #2 beta-globin; rename #3 myoglobin")
run(session, "compare #1,2,3 report globins_report")

run(session, "set bgColor white; lighting soft; graphics silhouettes true")
run(session, "hide atoms; show cartoons; view; turn y 30; view pad 0.02")
run(session, "save panel_deviation.png width 1200 height 1000 supersample 3")
run(session, "color #1 #8c8c8c; color #2 #1f77b4; color #3 #ff7f0e")
run(session, "save panel_species.png width 1200 height 1000 supersample 3")


def font(size, mono=False):
    names = ["consola.ttf", "cour.ttf"] if mono else ["arial.ttf", "DejaVuSans.ttf"]
    for n in names:
        try:
            return ImageFont.truetype(n, size)
        except OSError:
            pass
    return ImageFont.load_default()


# Alignment excerpt: first block of the multiple alignment
lines = open(os.path.join("globins_report", "multiple_alignment.txt"), encoding="utf-8").read().splitlines()
start = next(i for i, l in enumerate(lines) if l.startswith("alpha-globin "))
block = lines[start:start + 8]

W = 2400
def fit(path, box=(1200, 860)):
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, "white")
    from PIL import ImageChops
    bbox = ImageChops.difference(im, bg).getbbox()
    im = im.crop(bbox)
    scale = min(box[0] / im.width, box[1] / im.height) * 0.95
    im = im.resize((int(im.width * scale), int(im.height * scale)), Image.LANCZOS)
    canvas = Image.new("RGB", (box[0], 1000), "white")
    canvas.paste(im, ((box[0] - im.width) // 2, (box[1] - im.height) // 2 + 20))
    return canvas


species = fit("panel_species.png")
deviation = fit("panel_deviation.png")
plot = Image.open(os.path.join("globins_report", "deviation.png")).convert("RGB")
plot = plot.resize((W, int(plot.height * W / plot.width)))

mono = font(30, mono=True)
line_h = 38
aln_h = line_h * len(block) + 40
title_h = 70
H = title_h + 1000 + 30 + aln_h + plot.height
fig = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(fig)
head = font(34)
label = font(40)

d.text((20, 15), "a  Superposition, coloured by protein", fill="black", font=head)
d.text((W // 2 + 20, 15), "b  Coloured by Cα deviation from alpha-globin (blue = same, red = different)",
       fill="black", font=font(28))
fig.paste(species, (0, title_h))
fig.paste(deviation, (W // 2, title_h))
y = title_h + 1000 + 10
d.text((20, y), "c  Multiple alignment with secondary structure (H helix, E strand; * identical, = same SS in all)",
       fill="black", font=head)
y += 50
for l in block:
    d.text((40, y), l, fill="black", font=mono)
    y += line_h
y = title_h + 1000 + 30 + aln_h
fig.paste(plot, (0, y))
d.text((20, y + 5), "d", fill="black", font=label)

# Protein colour key on panel a
kx, ky = 40, title_h + 900
for name, col in (("alpha-globin (reference)", (140, 140, 140)), ("beta-globin", (31, 119, 180)),
                  ("myoglobin", (255, 127, 14))):
    d.rectangle([kx, ky + 8, kx + 30, ky + 38], fill=col)
    d.text((kx + 42, ky), name, fill="black", font=head)
    kx += 430

fig.save("compare_example.png", optimize=True)
for f in ("panel_species.png", "panel_deviation.png"):
    os.remove(f)
