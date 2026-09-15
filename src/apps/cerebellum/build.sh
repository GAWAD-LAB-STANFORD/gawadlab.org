#!/usr/bin/env bash
# Compile the cerebellum browser to WebAssembly and drop it at /cerebellum-atlas/.
#
# The data comes off Sherlock in three passes over the Carter et al. sqlite
# database (scripts in /oak/stanford/groups/cgawad/Scripts/atlas_export/):
#   export_cerebellum.py    per-cell t-SNE plus a 90x90 grid of mean log expression
#   export_cereb_groups.py  per-gene mean and detection rate in each cell type
#                           and each timepoint
#   reduce_cereb.py         re-bins to 42x42 and packs expression into a dense
#                           uint8 matrix - 25M sparse rows will not fit in a
#                           browser, 996 bins x 17,013 genes will
# prep_bundle.R then packs the lot into the single .rds the app reads.
#
#   bash src/apps/cerebellum/build.sh DATA_DIR [OUTDIR]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${1:?usage: build.sh DATA_DIR [OUTDIR]}"
OUT="${2:-$ROOT/cerebellum-atlas}"
STAGE="$(mktemp -d)/cerebellum-atlas"
mkdir -p "$STAGE"

cp "$SRC/app.R" "$STAGE/app.R"
Rscript "$SRC/prep_bundle.R" "$DATA" "$STAGE/bundle.rds"

Rscript -e "shinylive::export('$STAGE', '$OUT')"
echo "app.json: $(du -h "$OUT/app.json" | cut -f1)"
echo "built -> $OUT ($(du -sh "$OUT" | cut -f1))"
