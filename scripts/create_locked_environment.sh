#!/usr/bin/env bash
set -euo pipefail

component="${1:?Usage: create_locked_environment.sh main|idr|spp|explorer PREFIX}"
prefix="${2:?Provide a new environment prefix}"
case "$component" in
  main) name=omnomnomics ;;
  idr) name=omnomnomics-idr ;;
  spp) name=omnomnomics-spp ;;
  explorer) name=omnomnomics-explorer ;;
  *) printf 'Unknown environment: %s\n' "$component" >&2; exit 2 ;;
esac
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  printf 'These release locks support Linux x86-64 only.\n' >&2
  exit 2
fi
if [[ -e "$prefix" || -L "$prefix" ]]; then
  printf 'Refusing to modify an existing prefix: %s\n' "$prefix" >&2
  exit 2
fi
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
lock_dir="$repo/locks/linux-64"
test -s "$lock_dir/$name-lock.yml"
(cd "$lock_dir" && sha256sum --check --quiet SHA256SUMS)
micromamba_bin="${MICROMAMBA_EXE:-$(type -P micromamba || true)}"
test -x "$micromamba_bin"
"$micromamba_bin" --no-rc create --yes --prefix "$prefix" --file "$lock_dir/$name-lock.yml"
"$prefix/bin/python" "$repo/scripts/release/verify_environment.py" \
  --prefix "$prefix" --packages "$lock_dir/$name.packages.json"
