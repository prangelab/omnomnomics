#!/usr/bin/env bash
set -euo pipefail

MAIN_ENV="${OMNOMNOMICS_ENV:-omnomnomics}"
SPP_ENV="${OMNOMNOMICS_SPP_ENV:-omnomnomics-spp}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SPP_ENV_YAML="${REPO_ROOT}/environment.spp.yml"

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

if "${micromamba_bin}" env list | awk '{print $1}' | grep -qx "${SPP_ENV}"; then
    "${micromamba_bin}" install -y -n "${SPP_ENV}" -c conda-forge -c bioconda r-base=4.4 phantompeakqualtools=1.2.2
else
    "${micromamba_bin}" env create -y -f "${SPP_ENV_YAML}"
fi

main_prefix="$("${micromamba_bin}" run -n "${MAIN_ENV}" python -c "import sys; print(sys.prefix)")"
wrapper_path="${main_prefix}/bin/run_spp.R"

cat > "${wrapper_path}" <<EOF
#!/usr/bin/env bash
exec "${micromamba_bin}" run -n "${SPP_ENV}" run_spp.R "\$@"
EOF
chmod +x "${wrapper_path}"

"${micromamba_bin}" run -n "${SPP_ENV}" bash -lc "command -v run_spp.R >/dev/null"
"${micromamba_bin}" run -n "${MAIN_ENV}" bash -lc "command -v run_spp.R >/dev/null"
