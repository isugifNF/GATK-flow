#!/usr/bin/env bash
# CLI smoke checks for GATK-flow. Run from anywhere: bash tests/cli.sh
# Requires `nextflow` on PATH (CI pins NXF_VER=25.04.8).
set -uo pipefail

cd "$(dirname "$0")/.."

fail=0
pass() { echo "PASS: $1"; }
bad()  { echo "FAIL: $1"; fail=1; }

# Run a command, capturing combined output and exit code without aborting.
run() {
  rc=0
  out=$("$@" 2>&1) || rc=$?
}

# 1. Config parses and exposes manifest.version
run nextflow config -flat
if [[ $rc -eq 0 ]]; then pass "nextflow config"; else bad "nextflow config (exit $rc)"; echo "$out"; fi
version=$(sed -n "s/^manifest\.version = '\(.*\)'$/\1/p" <<<"$out")
if [[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.]+)?$ ]]; then
  pass "manifest.version is semver ($version)"
else
  bad "manifest.version is not semver: '$version'"
fi

# 2. --version
run nextflow run main.nf --version
lines=$(grep -c '^GATK-flow ' <<<"$out" || true)
[[ $rc -eq 0 ]] && pass "--version exits 0" || { bad "--version exit $rc"; echo "$out"; }
[[ $lines -eq 1 ]] && pass "--version prints one version line" || bad "--version printed $lines version lines"
grep -qx "GATK-flow $version" <<<"$out" && pass "--version matches manifest" || bad "--version does not print 'GATK-flow $version'"
grep -qiE 'genome|reads' <<<"$out" && bad "--version mentions missing genome/reads" || pass "--version needs no inputs"

# 3. --help
run nextflow run main.nf --help
[[ $rc -eq 0 ]] && pass "--help exits 0" || { bad "--help exit $rc"; echo "$out"; }
grep -q "GATK-flow $version" <<<"$out" && pass "--help shows version" || bad "--help missing 'GATK-flow $version'"
grep -q -- '--version' <<<"$out" && pass "--help documents --version" || bad "--help missing --version"

# 4. Unknown parameters are still rejected
run nextflow run main.nf --bogus 1
if [[ $rc -ne 0 ]] && grep -q 'not valid' <<<"$out"; then
  pass "unknown parameter rejected"
else
  bad "unknown parameter not rejected (exit $rc)"; echo "$out"
fi

exit $fail
