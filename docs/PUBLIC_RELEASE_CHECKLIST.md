# Public source release checklist

Before pushing to a public repository:

- [ ] Confirm supervisor, institutional and contributor approval for public code disclosure.
- [ ] Confirm the source code does not reveal findings you intend to embargo. Documentation is generic, but **the exact scripts may still contain project names, sample-count assertions, genomic coordinates and derived-result language**.
- [ ] Check for credentials, email addresses, local usernames, access tokens and private hostnames in all imported `.R` and `.py` files.
- [ ] Check for restricted data/metadata, result tables, figures, fitted models and manuscripts.
- [ ] Run `python3 tools/check_release_contents.py` and investigate all warnings.
- [ ] Confirm the list of tracked files with `git status` and `git ls-files` before pushing.
- [ ] Verify this is a **new repository without historic commits** that contained unpublished files.
- [ ] Choose and verify a software licence; do not imply open-source permission in its absence.
- [ ] If the repository becomes public, assume the committed source and history are permanently discoverable.
- [ ] Add an accurate software version and citation metadata once institutional ownership is confirmed.

Do not publish on the assumption that deleting a file later erases it from Git history or existing clones.
