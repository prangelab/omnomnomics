#!/usr/bin/env bash
set -euo pipefail

main="${1:?Usage: link_locked_helpers.sh MAIN_PREFIX IDR_PREFIX SPP_PREFIX}"
idr="${2:?Provide the IDR prefix}"
spp="${3:?Provide the SPP prefix}"
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
micromamba_bin="${MICROMAMBA_EXE:-$(type -P micromamba || true)}"
test -x "$micromamba_bin"
for component in main idr spp; do
  case "$component" in
    main) prefix="$main"; name=omnomnomics ;;
    idr) prefix="$idr"; name=omnomnomics-idr ;;
    spp) prefix="$spp"; name=omnomnomics-spp ;;
  esac
  "$main/bin/python" "$repo/scripts/release/verify_environment.py" \
    --prefix "$prefix" --packages "$repo/locks/linux-64/$name.packages.json"
done
for helper in idr run_spp.R; do
  if [[ "$helper" == idr ]]; then prefix="$idr"; else prefix="$spp"; fi
  wrapper="$(printf '#!/usr/bin/env bash\nexec %q run --prefix %q %q "$@"\n' "$micromamba_bin" "$prefix" "$helper")"
  if [[ -e "$main/bin/$helper" ]]; then
    if [[ "$(cat "$main/bin/$helper")" != "$wrapper" ]]; then
      printf 'Refusing to replace an existing executable: %s\n' "$main/bin/$helper" >&2
      exit 2
    fi
  else
    printf '%s\n' "$wrapper" > "$main/bin/$helper"
    chmod +x "$main/bin/$helper"
  fi
done
"$main/bin/idr" --version
"$micromamba_bin" run --prefix "$spp" Rscript --vanilla -e 'stopifnot(requireNamespace("spp", quietly=TRUE)); cat(as.character(packageVersion("spp")), "\n")'
