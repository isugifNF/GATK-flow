# GATK-flow Version Specification

## Goal

Give the pipeline a visible, maintainable release version so users can identify exactly which GATK-flow release they are running.

The user-facing version string will use this format:

```text
GATK-flow 0.1.0
```

The initial release version is `0.1.0`. The Nextflow manifest currently declares `1.0.0`, so implementation changes `manifest.version` to `0.1.0`. Later releases change that value as part of the release process.

## Current State

- `nextflow.config` currently declares `manifest.version = '1.0.0'`; it must be changed to `'0.1.0'`.
- `main.nf` supports `--help`, but it does not show the pipeline version.
- `main.nf` rejects parameters that are not in `parameters_valid`, so a new `--version` parameter must be explicitly allowed.
- Nextflow's own `-version` option reports the Nextflow runtime version, not the GATK-flow pipeline version. The pipeline option must therefore be `--version`.

## Version Contract

1. `manifest.version` in `nextflow.config` is the single source of truth.
2. Versions follow Semantic Versioning (`MAJOR.MINOR.PATCH`), for example `0.1.0`.
3. Running the following command prints the pipeline name and version, then exits successfully without requiring input files or starting processes:

   ```bash
   nextflow run main.nf --version
   ```

4. `--help` displays the same version and documents `--version`.
5. The version must not be copied into a second constant in `main.nf`; both outputs must read it from `workflow.manifest.version`.
6. `nextflow -version` and `nextflow run main.nf -version` remain Nextflow runtime behavior and are not used as the pipeline version interface.

## Implementation Steps

### 1. Define the CLI parameter default

In the `params` block in `nextflow.config`, add:

```groovy
version = false
```

Keep the existing manifest declaration as the canonical version:

```groovy
manifest {
  name = 'isugifNF/GATK'
  description = 'GATK pipeline'
  mainScript = 'main.nf'
  version = '0.1.0'
}
```

Do not introduce a separate version file or duplicate version constant.

### 2. Add a reusable version formatter

In `main.nf`, add a small function that obtains the value from the manifest and formats the public name consistently:

```groovy
def versionMsg() {
  "GATK-flow ${workflow.manifest.version}"
}
```

If Nextflow does not expose the configured manifest version in the target runtime (25.04.x per the project; 26.04.3 is installed locally and must also be verified), stop and verify the runtime API rather than adding an unrelated hard-coded fallback.

### 3. Show the version in help output

Add the formatted version near the beginning of `helpMsg()` so `--help` identifies the release being documented. Add this option to the arguments section:

```text
--version               Print the GATK-flow version and exit
```

The help output and standalone version output must use the same manifest value.

### 4. Implement the early-exit path

Before required input validation and before the `workflow` block can start work, handle the flag:

```groovy
if (params.version) {
  log.info versionMsg()
  exit 0
}
```

This check should be placed alongside the existing `params.help` early-exit handling. It must execute without `--genome`, reads, a profile, or a container runtime.

### 5. Allow the new parameter

Add `version` to `parameters_valid` in `main.nf`. This preserves the existing unknown-parameter validation while allowing `--version`.

### 6. Document version usage

Update the relevant README usage section, including the embedded `--help` output, with:

```bash
nextflow run main.nf --version
```

Explain that this reports the GATK-flow version, while `nextflow -version` reports the installed Nextflow runtime version.

### 7. Add focused checks

Add the version checks to the repository's test harness when it is created. At minimum, verify:

1. `nextflow config` succeeds.
2. `nextflow run main.nf --version` exits with status `0`.
3. Its output contains exactly one public version line matching `GATK-flow 0.1.0` for the current release.
4. The command does not complain about a missing genome or reads.
5. `nextflow run main.nf --help` exits with status `0` and contains both `GATK-flow 0.1.0` and `--version`.
6. An unknown long option is still rejected by the existing parameter validation.

Avoid asserting against unrelated Nextflow banner text, because that output may vary by runtime version.

## Release Procedure

For each release:

1. Choose the next Semantic Versioning value.
2. Change only `manifest.version` in `nextflow.config`.
3. Run the focused checks above and the existing smoke tests.
4. Confirm the README does not contain a stale hard-coded release value outside an intentional example.
5. Commit the version change with the release changes.
6. Create an annotated Git tag prefixed with `v`, matching the manifest version (for example, `v1.1.0`).
7. Push the commit and tag only after review.

Version increment rules:

- `PATCH`: backward-compatible bug fixes.
- `MINOR`: backward-compatible functionality.
- `MAJOR`: incompatible CLI, input, workflow, or output changes.

Pre-release versions may use Semantic Versioning suffixes such as `1.1.0-rc.1`; the corresponding Git tag is `v1.1.0-rc.1`.

## Files Expected to Change During Implementation

- `nextflow.config`: add the `version` parameter default and maintain `manifest.version`.
- `main.nf`: format, print, document, and allow `--version`.
- `README.md`: document how to query the pipeline version.
- Test harness file: add the focused CLI checks once the harness exists.

## Acceptance Criteria

The work is complete when all of the following are true:

- There is one canonical pipeline version value: `manifest.version`.
- `nextflow run main.nf --version` prints `GATK-flow <version>` and exits successfully.
- The version command does not require analysis inputs and starts no workflow processes.
- `--help` prints the same version and documents `--version`.
- Invalid parameters remain invalid.
- The documented release process keeps the manifest version and Git tag aligned.
