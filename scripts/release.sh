#!/bin/bash
set -e

# Default bump type
bump="patch"

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bump)
      bump="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Validate bump value
if [[ "$bump" != "patch" && "$bump" != "minor" && "$bump" != "major" ]]; then
  echo "Invalid bump type: $bump (must be patch, minor, or major)"
  exit 1
fi

# Git pull
git fetch origin main
git pull --rebase

# Ensure clean working tree
test -z "$(git status --porcelain)"

# Bump version
uv version --bump "$bump"

# Read version safely
version=$(grep -m1 version pyproject.toml | cut -d'"' -f2)

# Stage and commit
git add pyproject.toml uv.lock
git commit -m "Bump version to v$version"

# Tag release
git tag "v$version"

# Push commit + tag
git push origin main
git push origin "v$version"
