# Security and malware review, 2026-11-02

Scope: PyTransDecoder v0.2.1 (the PyPI sdist), Perl TransDecoder `master` (v7.0.0 plus README), the
TransDecoder wiki, and this paper repository and its manuscript submodule. All checks were read-only.
Tools are in `$P/env/secscan` (bandit, pip-audit, detect-secrets) and `$P/tools/gitleaks_8.21.2`
(release binary, checksum verified).

| Check | Result |
|---|---|
| gitleaks, full history of all branches of the 5 repos | no leaks |
| Credential-like filenames in history (.pypirc, .netrc, keys, gm_key) | none |
| GeneMarkS-T licence key content in any commit (pickaxe search, not printed) | none |
| bandit on PyTD 0.2.1 source | 9 low-severity findings, all `subprocess` use (below); no medium or high |
| pip-audit on the installed PyTD dependency set (biopython 1.88, numpy 2.4.6, tqdm 4.70.1) | no known vulnerabilities |
| PyTD: network, environment, home-directory, eval/exec, pickle, base64 | none (only documentation URLs and an SVG namespace) |
| PyTD bundled presets (35 files) | plain text and JSON; sequences are A/C/G/T/N only |
| Perl: network modules or download commands | none (URLs only in comments and an SVG namespace) |
| Perl: binary files in the tree | none apart from the gzipped presets |

External programs the code runs, all documented and user-requested:
- **PyTD:** gffread (genome mode); DIAMOND or makeblastdb/BLASTP, hmmpress/hmmsearch (homology options); samtools
  (genome mapping of BAM input). All are called with argument lists, never through a shell.
- **Perl:** LongOrfs/Predict steps, the decoy LongOrfs run, Pfam/BLAST helpers, and the `TransDecoder` wrapper,
  through `system()` with command strings built from user-supplied file names. This is the long-standing
  TransDecoder style, unchanged in kind since v6. A crafted file name could inject shell syntax, so it is a
  robustness issue, not malware.

Not done:
- An antivirus signature scan (ClamAV is not on conda-forge). It adds little for source code we wrote.
- Line-by-line human review of the changes since v6: about 1,900 Perl lines and the PyTD repository.

Recommended ongoing:
- **GitHub:** enable secret scanning with push protection, Dependabot alerts and CodeQL on both code repos.
- **PyPI:** 2FA on the account; a project-scoped token or trusted publishing from GitHub Actions; delete
  `gs://mdl-bhaas-tmp/.pypirc`.
- **Before each release:** rerun this checklist.
