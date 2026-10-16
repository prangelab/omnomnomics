#!/usr/bin/env bash
set -euo pipefail

MAIN_ENV="${OMNOMNOMICS_ENV:-omnomnomics}"
SPP_ENV="${OMNOMNOMICS_SPP_ENV:-omnomnomics-spp}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOCK_DIR="${REPO_ROOT}/locks/linux-64"

micromamba_bin="$(type -P micromamba || true)"
if [[ -z "${micromamba_bin}" ]]; then
    micromamba_bin="$(command -v micromamba || true)"
fi

if [[ -z "${micromamba_bin}" ]]; then
    echo "micromamba was not found on PATH." >&2
    exit 1
fi

if ! "${micromamba_bin}" run -n "${MAIN_ENV}" python -c "import sys" >/dev/null 2>&1; then
    echo "Main environment '${MAIN_ENV}' was not found. Create it from environment.yml first." >&2
    exit 1
fi

main_prefix="$("${micromamba_bin}" run -n "${MAIN_ENV}" python -c "import sys; print(sys.prefix)")"
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "The release SPP lock supports Linux x86-64 only." >&2
    exit 1
fi
(cd "${LOCK_DIR}" && sha256sum --check --quiet SHA256SUMS)
if ! "${micromamba_bin}" run -n "${SPP_ENV}" python -c 'import sys' >/dev/null 2>&1; then
    "${micromamba_bin}" --no-rc create -y -n "${SPP_ENV}" -f "${LOCK_DIR}/omnomnomics-spp-lock.yml"
fi
helper_prefix="$("${micromamba_bin}" run -n "${SPP_ENV}" python -c 'import sys; print(sys.prefix)')"
"${main_prefix}/bin/python" "${REPO_ROOT}/scripts/release/verify_environment.py" \
    --prefix "${helper_prefix}" --packages "${LOCK_DIR}/omnomnomics-spp.packages.json"
wrapper_path="${main_prefix}/bin/run_spp.R"

cat > "${wrapper_path}" <<EOF
#!/usr/bin/env bash
exec "${micromamba_bin}" run -n "${SPP_ENV}" run_spp.R "\$@"
EOF
chmod +x "${wrapper_path}"

"${micromamba_bin}" run -n "${SPP_ENV}" Rscript --vanilla -e 'stopifnot(requireNamespace("spp", quietly=TRUE))'
"${micromamba_bin}" run -n "${MAIN_ENV}" python -c 'import shutil; assert shutil.which("run_spp.R")'
