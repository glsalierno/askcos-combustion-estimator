#!/usr/bin/env bash
# Log process + latest iterative CSV stats every INTERVAL seconds (default 600 = 10 min).
# Env: PROJECT_ROOT, OUTPUT_GLOB, INPUT_CSV, CMD_PATTERN, MONITOR_LOG.
set -u

INTERVAL="${1:-600}"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
LOG="${MONITOR_LOG:-$PROJECT_ROOT/monitor_run.log}"
OUTPUT_GLOB="${OUTPUT_GLOB:-$PROJECT_ROOT/output_iterative_n100_run_*.csv}"
INPUT_CSV="${INPUT_CSV:-$PROJECT_ROOT/batch_smiles_n100.csv}"
CMD_PATTERN="${CMD_PATTERN:-process_reactivity_oxygen.py}"

log() {
  local line="[$(date -Iseconds)] $*"
  printf '%s\n' "$line" >>"$LOG"
  if [[ -t 1 ]]; then printf '%s\n' "$line"; fi
}

total_inputs=100
if [[ -f "$INPUT_CSV" ]]; then
  total_inputs=$(( $(wc -l <"$INPUT_CSV") - 1 ))
fi

pick_latest_csv() {
  local latest="" m=0 t f
  shopt -s nullglob
  for f in $OUTPUT_GLOB; do
    [[ -f "$f" ]] || continue
    if stat --version >/dev/null 2>&1; then
      t=$(stat -c %Y "$f")
    else
      t=$(stat -f %m "$f")
    fi
    if (( t > m )); then m=$t; latest=$f; fi
  done
  shopt -u nullglob
  printf '%s' "$latest"
}

while true; do
  log "=== monitor tick ==="
  if out=$(pgrep -af "$CMD_PATTERN" 2>/dev/null); then
    n=$(printf '%s\n' "$out" | wc -l)
    log "process: RUNNING ($n lines from pgrep -af)"
    while IFS= read -r line; do log "  $line"; done <<<"$out"
  else
    log "process: NOT RUNNING"
  fi

  latest="$(pick_latest_csv)"
  if [[ -n "$latest" ]]; then
    log "latest_csv: $latest"
    if stat --version >/dev/null 2>&1; then
      log "  $(stat -c 'mtime=%y bytes=%s' "$latest")"
    else
      log "  $(stat -f 'mtime=%Sm bytes=%z' "$latest")"
    fi
    stats=$(MONITOR_CSV="$latest" MONITOR_N="$total_inputs" python3 - <<'PY'
import os
import pandas as pd

p = os.environ["MONITOR_CSV"]
n_in = int(os.environ["MONITOR_N"])
df = pd.read_csv(p)
n = df["original_smiles"].nunique()
mx = float(df["probability"].max())
md = float(df["cumulative_delta_g_kj_per_mol"].min())
print(
    f"rows={len(df)} distinct_inputs={n}/{n_in} "
    f"max_prob={mx:.6g} min_delta_g_kj={md:.6g}"
)
PY
)
    log "  $stats"
  else
    log "latest_csv: (none matching OUTPUT_GLOB)"
  fi

  sleep "$INTERVAL"
done
