#!/usr/bin/env bash
set -euo pipefail

MAIN_ENV="${OMNOMNOMICS_ENV:-omnomnomics}"
IDR_ENV="${OMNOMNOMICS_IDR_ENV:-omnomnomics-idr}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IDR_ENV_YAML="${REPO_ROOT}/environment.idr.yml"

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

if "${micromamba_bin}" env list | awk '{print $1}' | grep -qx "${IDR_ENV}"; then
    "${micromamba_bin}" install -y -n "${IDR_ENV}" -c conda-forge -c bioconda python=3.10 idr "numpy<1.24"
else
    "${micromamba_bin}" env create -f "${IDR_ENV_YAML}"
fi

main_prefix="$("${micromamba_bin}" run -n "${MAIN_ENV}" python -c "import sys; print(sys.prefix)")"
wrapper_path="${main_prefix}/bin/idr"

cat > "${wrapper_path}" <<EOF
#!/usr/bin/env bash
exec "${micromamba_bin}" run -n "${IDR_ENV}" idr "\$@"
EOF
chmod +x "${wrapper_path}"

"${micromamba_bin}" run -n "${MAIN_ENV}" idr --version
