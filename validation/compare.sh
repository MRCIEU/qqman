#!/usr/bin/env bash
# Renders the cases in reference.R with qqman-rs, as PDF and PNG, and compares
# each against R/qqman's output pixel by pixel.
#
#   validation/compare.sh                 # compare against existing reference/
#   validation/compare.sh --reference     # rerun reference.R first (needs R + qqman)
#
# PDFs: both are rasterised with Ghostscript at DPI (default 100). A pixel
# differs when any channel differs by more than THRESHOLD/255 (default 60).
#
# PNGs (150 dpi) are compared directly, with a one-pixel
# tolerance: a pixel counts only if no pixel within one pixel of it in the
# other image matches it (within THRESHOLD), checked both ways. R's PNG device
# anti-aliases and draws glyphs differently from rs-graphics (and on macOS,
# quartz also colour-manages: gray60 comes out darker), so a plain count is
# dominated by edge noise. The tolerant count still catches a missing or moved
# element. On 19 Oct 2026, at the default sizes, the baselines were 556-1350
# pixels, while removing the threshold lines gave 8655, dropping the
# annotation labels 4592, and points at cex 0.8 gave 5033. PDF baselines were
# 29-133, and the same two regressions gave 5209 and 1585. reference/png-device.txt records which device made the reference
# PNGs; references made with cairo will differ.
#
# Exits 1 if any case has more differing pixels than PDF_MAX_DIFF (default
# 200) or PNG_MAX_DIFF (default 2000).
#
# Needs: cargo, gs, python3.
set -euo pipefail

DPI=${DPI:-100}
THRESHOLD=${THRESHOLD:-60}
PDF_MAX_DIFF=${PDF_MAX_DIFF:-200}
PNG_MAX_DIFF=${PNG_MAX_DIFF:-2000}

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$(dirname "$here")
data=$here/data
ref=$here/reference
out=$here/out
cases=(manhattan manhattan_main manhattan_chr3 manhattan_highlight manhattan_annotate qq)

for tool in cargo gs python3; do
    command -v "$tool" >/dev/null || { echo "compare.sh: $tool not found" >&2; exit 2; }
done

if [[ ${1:-} == --reference ]]; then
    Rscript "$here/reference.R"
fi
for name in "${cases[@]}"; do
    for ext in pdf png; do
        if [[ ! -f $ref/$name.$ext ]]; then
            echo "compare.sh: $ref/$name.$ext missing; run with --reference (needs R + qqman)" >&2
            exit 2
        fi
    done
done

(cd "$root" && cargo build --release --quiet)
qqman=$root/target/release/qqman

# The same cases as reference.R, in each format.
mkdir -p "$out"
awk -F'\t' 'NR == 1 || $2 == 3' "$data/gwasResults.tsv" > "$out/chr3.tsv"
# Rendered at the CLI's default sizes (manhattan 14 x 7 in, qq 7 x 7 in, PNG
# at 150 dpi), which reference.R matches. If the defaults change, the size check
# below fails until reference.R is updated and rerun.
render() {
    local name=$1; shift
    for ext in pdf png; do
        "$qqman" "$@" -o "$out/$name.$ext" 2>/dev/null
    done
}
render manhattan           manhattan "$data/gwasResults.tsv"
render manhattan_main      manhattan "$data/gwasResults.tsv" --main "Manhattan Plot" --ylim 0,10 \
                           --cex 0.6 --cex-axis 0.9 --col blue4,orange3 \
                           --suggestiveline false --genomewideline false
render manhattan_chr3      manhattan "$out/chr3.tsv"
render manhattan_highlight manhattan "$data/gwasResults.tsv" --highlight-file "$data/snpsOfInterest.txt"
render manhattan_annotate  manhattan "$data/gwasResults.tsv" --annotate-pval 0.01
render qq                  qq "$data/gwasResults.tsv"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
pairs=()
for name in "${cases[@]}"; do
    gs -q -sDEVICE=ppmraw -r"$DPI" -o "$tmp/r_$name.ppm" "$ref/$name.pdf"
    gs -q -sDEVICE=ppmraw -r"$DPI" -o "$tmp/rs_$name.ppm" "$out/$name.pdf"
    pairs+=("$name.pdf" "$tmp/r_$name.ppm" "$tmp/rs_$name.ppm")
done
for name in "${cases[@]}"; do
    pairs+=("$name.png" "$ref/$name.png" "$out/$name.png")
done

python3 -I - "$THRESHOLD" "$PDF_MAX_DIFF" "$PNG_MAX_DIFF" "${pairs[@]}" <<'EOF'
import struct
import sys
import zlib


def load_ppm(path):
    """(width, height, RGB bytes) of a binary P6 PPM, skipping header comments."""
    data = open(path, "rb").read()
    fields, pos = [], 2  # after "P6"
    while len(fields) < 3:
        while data[pos:pos + 1].isspace():
            pos += 1
        if data[pos:pos + 1] == b"#":
            pos = data.index(b"\n", pos) + 1
            continue
        end = pos
        while not data[end:end + 1].isspace():
            end += 1
        fields.append(int(data[pos:end]))
        pos = end
    w, h, _maxval = fields
    return w, h, data[pos + 1:pos + 1 + w * h * 3]


def load_png(path):
    """(width, height, RGB bytes) of an 8-bit, non-interlaced RGB or RGBA PNG,
    with any alpha composited onto white."""
    data = open(path, "rb").read()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"{path}: not a PNG")
    pos, idat = 8, []
    while pos < len(data):
        length, kind = struct.unpack(">I4s", data[pos:pos + 8])
        body = data[pos + 8:pos + 8 + length]
        if kind == b"IHDR":
            w, h, depth, colour, _, _, interlace = struct.unpack(">IIBBBBB", body)
            if depth != 8 or colour not in (2, 6) or interlace:
                raise ValueError(f"{path}: unsupported PNG (depth {depth}, colour {colour}, interlace {interlace})")
            bpp = 3 if colour == 2 else 4
        elif kind == b"IDAT":
            idat.append(body)
        elif kind == b"IEND":
            break
        pos += 12 + length
    raw = zlib.decompress(b"".join(idat))
    stride = w * bpp
    prev = bytearray(stride)
    rgb = bytearray(w * h * 3)
    for y in range(h):
        start = y * (stride + 1)
        kind, line = raw[start], bytearray(raw[start + 1:start + 1 + stride])
        if kind == 1:  # Sub
            for i in range(bpp, stride):
                line[i] = (line[i] + line[i - bpp]) & 255
        elif kind == 2:  # Up
            for i in range(stride):
                line[i] = (line[i] + prev[i]) & 255
        elif kind == 3:  # Average
            for i in range(stride):
                left = line[i - bpp] if i >= bpp else 0
                line[i] = (line[i] + ((left + prev[i]) >> 1)) & 255
        elif kind == 4:  # Paeth
            for i in range(stride):
                a = line[i - bpp] if i >= bpp else 0
                b = prev[i]
                c = prev[i - bpp] if i >= bpp else 0
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                line[i] = (line[i] + (a if pa <= pb and pa <= pc else b if pb <= pc else c)) & 255
        out = y * w * 3
        if bpp == 3:
            rgb[out:out + w * 3] = line
        else:
            for x in range(w):
                r, g, b, alpha = line[4 * x:4 * x + 4]
                for k, v in enumerate((r, g, b)):
                    rgb[out + 3 * x + k] = (v * alpha + 255 * (255 - alpha) + 127) // 255
        prev = line
    return w, h, bytes(rgb)


def plain_count(a, b, w, h, t):
    """Pixels whose colour differs by more than t in any channel."""
    return sum(
        1 for i in range(0, len(a), 3)
        if max(abs(a[i] - b[i]), abs(a[i + 1] - b[i + 1]), abs(a[i + 2] - b[i + 2])) > t
    )


def tolerant_count(a, b, w, h, t):
    """Pixels of a with no pixel of b within one pixel (8 neighbours or
    itself) within t, plus the same the other way round."""
    def one_way(a, b):
        n = 0
        for p in range(w * h):
            i = 3 * p
            if max(abs(a[i] - b[i]), abs(a[i + 1] - b[i + 1]), abs(a[i + 2] - b[i + 2])) <= t:
                continue
            y, x = divmod(p, w)
            if not any(
                max(abs(a[i] - b[j]), abs(a[i + 1] - b[j + 1]), abs(a[i + 2] - b[j + 2])) <= t
                for yy in range(max(0, y - 1), min(h, y + 2))
                for xx in range(max(0, x - 1), min(w, x + 2))
                for j in (3 * (yy * w + xx),)
            ):
                n += 1
        return n
    return one_way(a, b) + one_way(b, a)


threshold, pdf_max, png_max = (int(v) for v in sys.argv[1:4])
args = sys.argv[4:]
failed = False
print(f"{'case':<26}{'differing px':>14}{'of':>10}{'limit':>8}")
for name, r_path, rs_path in zip(args[0::3], args[1::3], args[2::3]):
    is_png = name.endswith(".png")
    load = load_png if is_png else load_ppm
    limit = png_max if is_png else pdf_max
    w, h, a = load(r_path)
    w2, h2, b = load(rs_path)
    if (w, h) != (w2, h2):
        print(f"{name:<26}  size differs: R {w}x{h}, qqman-rs {w2}x{h2}")
        failed = True
        continue
    n = (tolerant_count if is_png else plain_count)(a, b, w, h, threshold)
    flag = "  <-- over limit" if n > limit else ""
    failed |= n > limit
    print(f"{name:<26}{n:>14}{w * h:>10}{limit:>8}{flag}")
print(f"(pixels differing by more than {threshold}/255; PNGs with a one-pixel tolerance)")
sys.exit(1 if failed else 0)
EOF
