# Security Policy

## Supported Versions

plotlyMol is pre-1.0 software (currently `0.2.x`) under active development. There is no
long-term-support branch: **only the latest release published on PyPI** receives security
fixes.

| Version        | Supported          |
| -------------- | ------------------ |
| latest (0.2.x) | :white_check_mark: |
| < 0.2          | :x:                |

If you're on an older version, please upgrade before reporting an issue that a newer
release may have already fixed.

## Reporting a Vulnerability

Please **do not open a public GitHub issue** for security vulnerabilities.

Instead, report it privately using GitHub's built-in reporting tool:

1. Go to the repository's [Security tab](../../security/advisories).
2. Click **"Report a vulnerability"** to open a private advisory.

This reaches the maintainers directly without disclosing the issue publicly before a fix
is available. If you're unable to use GitHub's private reporting, open an issue that omits
exploit details and asks a maintainer to follow up privately instead.

When reporting, please include:

- A description of the vulnerability and its potential impact.
- Steps to reproduce it (a minimal script or molecule input is ideal).
- The plotlyMol version and Python version you're using.

### What to expect

- We'll acknowledge new reports within about a week.
- If accepted, we'll work on a fix and coordinate a release with you, crediting you in the
  release notes unless you'd rather stay anonymous.
- If declined (not reproducible, out of scope, etc.), we'll explain why.

### Scope

This covers the plotlyMol source in this repository. Vulnerabilities in dependencies
(RDKit, Plotly, NumPy, Kaleido, etc.) should generally go to those projects directly —
happy to help route a report if you're not sure where it belongs.
