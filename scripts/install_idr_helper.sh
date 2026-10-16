#!/usr/bin/env bash
set -euo pipefail

MAIN_ENV="${OMNOMNOMICS_ENV:-omnomnomics}"
IDR_ENV="${OMNOMNOMICS_IDR_ENV:-omnomnomics-idr}"

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
    echo "The release IDR lock supports Linux x86-64 only." >&2
    exit 1
fi
(cd "${LOCK_DIR}" && sha256sum --check --quiet SHA256SUMS)
if ! "${micromamba_bin}" run -n "${IDR_ENV}" python -c 'import sys' >/dev/null 2>&1; then
    "${micromamba_bin}" --no-rc create -y -n "${IDR_ENV}" -f "${LOCK_DIR}/omnomnomics-idr-lock.yml"
fi
helper_prefix="$("${micromamba_bin}" run -n "${IDR_ENV}" python -c 'import sys; print(sys.prefix)')"
"${main_prefix}/bin/python" "${REPO_ROOT}/scripts/release/verify_environment.py" \
    --prefix "${helper_prefix}" --packages "${LOCK_DIR}/omnomnomics-idr.packages.json"
wrapper_path="${main_prefix}/bin/idr"

cat > "${wrapper_path}" <<EOF
#!/usr/bin/env bash
exec "${micromamba_bin}" run -n "${IDR_ENV}" idr "\$@"
EOF
chmod +x "${wrapper_path}"

"${micromamba_bin}" run -n "${MAIN_ENV}" idr --version
