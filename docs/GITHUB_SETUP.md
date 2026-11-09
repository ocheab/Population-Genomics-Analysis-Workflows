# Publish code with generic documentation

Create a new, empty GitHub repository. Choose **private** if you have not yet cleared public code release with your supervisor/institution. Do not initialise the remote repository with an extra README to avoid a first-commit conflict.

From this clean source directory:

```bash
# Optional: import unmodified Python files from your local project first.
bash tools/import_python_from_local_project.sh /absolute/path/to/project/scripts
python3 tools/check_release_contents.py

# Review the paths before first commit.
git init -b main
git status --short
# Only add source and approved documentation:
git add README.md .gitignore .gitattributes R docs tools python
# Review staged contents:
git diff --cached --stat
git diff --cached --name-only
git commit -m "Publish analysis source and generic methodology documentation"

git remote add origin https://github.com/USERNAME/REPOSITORY.git
git push -u origin main
```

Replace `USERNAME` and `REPOSITORY` with the actual GitHub account and new repository name. The example URL is a placeholder, not an existing repository. Configure Git credentials through your preferred authentication method before pushing. The files excluded by `.gitignore` remain excluded even if later generated locally, but `.gitignore` does **not** remove files that are already tracked by Git.

For later updates, use tagged releases after review. If you need to keep findings unpublished, leave the repository private or publish only the source that you explicitly approve.
