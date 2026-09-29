#!/usr/bin/env bash
set -euo pipefail
workflow_dir="$(cd "$(dirname "$0")" && pwd)"
step=${1:?选择分析步骤：10_paga|11_scvelo}
shift
case "$step" in
  10_paga|11_scvelo) python "$workflow_dir/$step.py" "$@" ;;
  *) printf '%s\n' '未知步骤' >&2; exit 2 ;;
esac
