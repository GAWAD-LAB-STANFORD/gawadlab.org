#!/usr/bin/env bash
# Compile the pediatric ALL treatment-resistance browser to WebAssembly and drop
# it at /pall-resistance/.
#
#   bash src/apps/pall/build.sh DATA_DIR [OUTDIR]
#
# DATA_DIR is the manuscript's 08_Source_Data tree. prep_bundle.R reads the
# figure-level tables out of it and writes the single .rds the app loads.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${1:?usage: build.sh DATA_DIR [OUTDIR]}"
OUT="${2:-$ROOT/pall-resistance}"
STAGE="$(mktemp -d)/pall-resistance"
mkdir -p "$STAGE"

cp "$SRC/app.R" "$STAGE/app.R"
Rscript "$SRC/prep_bundle.R" "$DATA" "$STAGE/bundle.rds"

Rscript -e "shinylive::export('$STAGE', '$OUT')"
echo "app.json: $(du -h "$OUT/app.json" | cut -f1)"
echo "built -> $OUT ($(du -sh "$OUT" | cut -f1))"
