#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 || ! -d "$1" ]]; then
  echo "Usage: bash tools/import_python_from_local_project.sh /absolute/path/to/original/project/scripts" >&2
  exit 2
fi
source_dir=$(cd "$1" && pwd -P)
repo_dir=$(cd "$(dirname "$0")/.." && pwd -P)
dest_dir="$repo_dir/python"
if [[ "$source_dir" == "$dest_dir" ]]; then
  echo "Source and destination are identical" >&2; exit 2
fi
mkdir -p "$dest_dir"
shopt -s nullglob
files=("$source_dir"/*.py)
if (( ${#files[@]} == 0 )); then
  echo "No .py files in $source_dir" >&2; exit 1
fi
for f in "${files[@]}"; do
  cp -- "$f" "$dest_dir/$(basename "$f")"
  printf 'IMPORTED (byte-identical): %s\n' "$(basename "$f")"
done
(
  cd "$repo_dir"
  for f in python/*.py; do sha256sum "$f"; done
)
echo "Review all imported scripts for secrets and unpublished findings before git push."
