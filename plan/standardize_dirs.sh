#!/usr/bin/env bash
# standardize_dirs.sh
# Standardizes folder structure across all research experiment directories.
# Usage:
#   ./standardize_dirs.sh --dry-run   # preview actions without changing anything
#   ./standardize_dirs.sh             # apply changes

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DRY_RUN=false

if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "=== DRY RUN — no files will be changed ==="
fi

# ── helpers ──────────────────────────────────────────────────────────────────

do_mkdir() {
    local dir="$1"
    if [[ "$DRY_RUN" == true ]]; then
        echo "  [mkdir] $dir"
    else
        mkdir -p "$dir"
        touch "$dir/.gitkeep"
    fi
}

do_move() {
    local src="$1" dst="$2"
    if [[ "$DRY_RUN" == true ]]; then
        echo "  [mv]    $src  →  $dst"
    else
        mv "$src" "$dst"
    fi
}

do_write_license() {
    local path="$1"
    if [[ "$DRY_RUN" == true ]]; then
        echo "  [LICENSE] $path"
    else
        cat > "$path" << 'EOF'
MIT License

Copyright (c) 2024 Gian Nico Benucci

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOF
    fi
}

# ── config ────────────────────────────────────────────────────────────────────

# These two dirs keep Suppl_File_X flat naming — create skeleton only, no file moves
SUPPL_DIRS=(
    "Benucci_etal_2020_Mao-Tofu_Microbiome"
    "Benucci_etal_2020_Patient_Propagules"
)

# These dirs already have internal structure — create missing folders + LICENSE only
PARTIAL_DIRS=(
    "Benucci_etal_2018_FijiSoilMicrobiome"
    "Benucci_etal_2019_Lycopods"
    "Benucci_etal_2019_MorchellaMicrobiome"
    "Benucci_etal_2025_CompetitionReleaseTruffle"
    "Garcia-Barreda_et_al_TruffleNestBacteria"
    "Kelley_etal_2024_ArabidopsisDrought"
    "Liu_etal_2026_AMF_FertUnfertSwitchgrass"
)

# File sorting rules — extension → destination folder (relative to experiment dir)
# Order matters: more-specific patterns first
declare -A EXT_MAP=(
    [R]="code"
    [r]="code"
    [Rmd]="code"
    [rmd]="code"
    [sh]="code"
    [rds]="datasets"
    [RDS]="datasets"
    [RData]="datasets"
    [Rdata]="datasets"
    [rda]="datasets"
    [csv]="datasets"
    [tsv]="datasets"
    [txt]="datasets"
    [xlsx]="datasets"
    [xls]="datasets"
    [biom]="datasets"
    [fasta]="datasets"
    [fa]="datasets"
    [fastq]="datasets"
    [pdf]="misc"
    [eps]="misc"
    [png]="misc"
    [jpg]="misc"
    [jpeg]="misc"
)

# Files that must stay at the experiment root, never moved
ROOT_FILES=("README.md" "LICENSE" "LICENSE.txt" "*.Rproj" "renv.lock")

# ── functions ─────────────────────────────────────────────────────────────────

is_root_file() {
    local f="$1"
    local base
    base="$(basename "$f")"
    [[ "$base" == README.md ]] && return 0
    [[ "$base" == LICENSE* ]] && return 0
    [[ "$base" == *.Rproj ]] && return 0
    [[ "$base" == renv.lock ]] && return 0
    return 1
}

ensure_skeleton() {
    local expdir="$1"
    echo "  skeleton: code/ datasets/ functions/ misc/"
    for folder in code datasets functions misc; do
        [[ -d "$expdir/$folder" ]] || do_mkdir "$expdir/$folder"
    done
}

ensure_license() {
    local expdir="$1"
    if [[ ! -f "$expdir/LICENSE" && ! -f "$expdir/LICENSE.txt" ]]; then
        echo "  license: adding MIT LICENSE"
        do_write_license "$expdir/LICENSE"
    fi
}

sort_loose_files() {
    local expdir="$1"
    # Only process regular files directly at the experiment root
    while IFS= read -r -d '' f; do
        local base ext dst
        base="$(basename "$f")"
        is_root_file "$f" && continue
        ext="${base##*.}"
        dst="${EXT_MAP[$ext]:-}"
        if [[ -n "$dst" ]]; then
            do_move "$f" "$expdir/$dst/$base"
        else
            echo "  [skip]  $base (no rule for extension .$ext)"
        fi
    done < <(find "$expdir" -maxdepth 1 -type f -print0)
}

# ── main loop ─────────────────────────────────────────────────────────────────

for expdir in "$REPO_ROOT"/*/; do
    expname="$(basename "$expdir")"

    # Skip hidden dirs and plan/ itself
    [[ "$expname" == .* ]] && continue
    [[ "$expname" == "plan" ]] && continue

    echo ""
    echo "━━━ $expname"

    # Determine which category this dir falls into
    is_suppl=false
    is_partial=false

    for d in "${SUPPL_DIRS[@]}"; do
        [[ "$expname" == "$d" ]] && { is_suppl=true; break; }
    done

    for d in "${PARTIAL_DIRS[@]}"; do
        [[ "$expname" == "$d" ]] && { is_partial=true; break; }
    done

    ensure_skeleton "$expdir"
    ensure_license  "$expdir"

    if [[ "$is_suppl" == true ]]; then
        echo "  [Suppl_File dir] — skeleton created, files left in place"
        continue
    fi

    if [[ "$is_partial" == true ]]; then
        echo "  [partial dir] — folders ensured, no file moves"
        continue
    fi

    # Flat dir: sort loose files
    sort_loose_files "$expdir"
done

echo ""
echo "=== Done ==="
