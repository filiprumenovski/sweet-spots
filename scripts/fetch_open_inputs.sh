#!/usr/bin/env bash

set -euo pipefail

repository_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
archive_name="sweet-spots-open-inputs-v1.tar.gz"
archive_url="https://github.com/filiprumenovski/sweet-spots/releases/download/inputs-v1/${archive_name}"
expected_sha256="4010a15610f2a8d62c5f6ae0728249bf6f21353e42477ebfd359cb8d4c0fd3e8"

command -v curl >/dev/null || {
    printf '%s\n' "curl is required to download ${archive_name}." >&2
    exit 1
}
command -v sha256sum >/dev/null || {
    printf '%s\n' "sha256sum is required to verify ${archive_name}." >&2
    exit 1
}

archive_path="$(mktemp "${TMPDIR:-/tmp}/sweet-spots-inputs.XXXXXX.tar.gz")"
cleanup() {
    rm -f -- "${archive_path}"
}
trap cleanup EXIT

printf 'Downloading %s\n' "${archive_url}"
curl --fail --location --retry 3 --output "${archive_path}" "${archive_url}"
printf '%s  %s\n' "${expected_sha256}" "${archive_path}" | sha256sum --check

printf 'Extracting into %s\n' "${repository_root}/inputs"
tar -xzf "${archive_path}" -C "${repository_root}"

missing=0
for relative_path in \
    analysis/revalidation/data/atlas_unambiguous.csv \
    data/external/multispecies_oglcnac/rice_sequences.fasta \
    data/processed/landscape/ptm_unified.parquet
do
    if [[ ! -f "${repository_root}/inputs/${relative_path}" ]]; then
        printf 'User-supplied input still required: inputs/%s\n' "${relative_path}"
        missing=1
    fi
done

if [[ "${missing}" -eq 0 ]]; then
    printf '%s\n' "The complete frozen input tree is present."
else
    printf '%s\n' "Open inputs installed. See docs/input_data.md for the three licensed inputs."
fi
