# Releasing SVscanner

The release number lives in exactly one place: the [`VERSION`](../VERSION) file at the
repository root, as a bare `X.Y.Z`. Edit it by hand when cutting a release.

| Where the version appears | How it stays correct |
|---|---|
| `VERSION` | Source of truth, edited by hand. |
| `scripts/run_workflow.sh` (`--version`, run banner) | Reads `VERSION` at startup. Never hardcode it. |
| git tag `vX.Y.Z` | Must match `VERSION` at the tagged commit — checked by CI and by the pre-push hook. |
| Container image tag on GHCR | Built and pushed from the tag by [`docker-publish.yml`](../.github/workflows/docker-publish.yml). Nothing to edit. |
| `README.md`, `docs/docker.md` (`ghcr.io/…/svscanner:X.Y.Z` examples) | Pinned versions in the usage examples. Checked as a warning only. |

## Cutting a release

**Bump `VERSION` before the release commit is merged** — the tag is only correct if the
commit it points at already carries the new number.

```bash
# 1. Edit VERSION to the new number and commit it.
vi VERSION
git commit -am "release v0.6.1"
git push origin dev
```

2. Merge `dev` into `main` as usual.

3. Create the tag and release on the GitHub website
   (*Releases → Draft a new release → Choose a tag → Create new tag*), targeting that
   merge commit. Name the tag `v0.6.1`.

4. **Check the `version check` run** under the Actions tab. Publishing the release
   fires it against the new tag; if it fails, the tag disagrees with `VERSION` and the
   release needs re-cutting (delete the release and its tag, fix `VERSION`, tag again).

5. **Check the `docker publish` run** under the Actions tab. It builds the image, runs
   its smoke tests (including RepeatMasker on an SVA in both orientations) and pushes
   `ghcr.io/gentechgp/svscanner:X.Y.Z`, `:X.Y` and `:latest`. The first push of a package
   creates it as **private** — set it to public once, under
   *Packages → svscanner → Package settings*, or nobody else can pull it.

6. Pull the published image and confirm it reports what you expect:

```bash
singularity pull svscanner_0.6.1.sif docker://ghcr.io/gentechgp/svscanner:0.6.1
singularity exec svscanner_0.6.1.sif svscanner --version   # must print: SVscanner v0.6.1
```

7. If the README or [docs/docker.md](docker.md) pin a `ghcr.io/…/svscanner:X.Y.Z`
   example, update it to the new version.

## The two checks

Both run [`check_version.sh`](../scripts/check_version.sh), which reads files straight
out of a commit — nothing is checked out and the working tree is untouched.

**[CI](../.github/workflows/version-check.yml)** — runs on every branch push, pull
request and tag push. This is the one that covers tags created on the GitHub website,
which never reach a local hook. It reports *after* the tag exists, so treat a red
`version check` on a tag as "re-cut this release".

**[`.githooks/pre-push`](../.githooks/pre-push)** — refuses the push before anything
leaves your machine, for tags pushed from the command line. Enable once per clone:

```bash
git config core.hooksPath .githooks
```

`core.hooksPath` is local to each clone, so a fresh clone needs this again.
`git push --no-verify` bypasses it.

### What counts as a failure

Hard errors, which fail CI and block a local push:

- a tag `vX.Y.Z` published on a commit whose `VERSION` does not say `X.Y.Z` — the
  mistake that produced a v0.6.0 tag on a 0.5.2 script;
- `VERSION` missing or malformed;
- `run_workflow.sh` hardcoding a version instead of reading `VERSION`.

Warning only, never blocking:

- a `ghcr.io/…/svscanner:X.Y.Z` example in the README or `docs/docker.md` differing from
  `VERSION`. An image *is* published for every tag, so this one is genuinely stale rather
  than legitimately lagging — but the check necessarily runs before the image it names
  exists, so blocking a release on it would be circular.

Commits from before the `VERSION` file existed are skipped, so old branches still push.

## Checking by hand

```bash
scripts/check_version.sh                 # the working tree
scripts/check_version.sh --rev <commit>  # any commit, without checking it out
scripts/check_version.sh --expect 0.6.1  # as CI does for tag v0.6.1
```

Only `ghcr.io/…/svscanner:X.Y.Z` image references are compared, so prose that names an
older release, or the `X.Y.Z` placeholder in the usage examples, is left alone.
