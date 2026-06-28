#!/usr/bin/env bash
set -euo pipefail

echo "[1/4] System dependencies..."
sudo apt-get update -y
sudo apt-get install -y \
  curl \
  ca-certificates \
  zstd \
  build-essential \
  libgl1

echo "[2/4] Install uv..."
curl -Ls https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.cargo/bin:$PATH"

echo "[3/4] Create and sync Python environment..."
# Create virtual environment and install from pyproject.toml file
uv sync --group dev
# Activate virtual environment
source .venv/bin/activate

echo "[4/4] Setup Jupyter kernel..."
python -m ipykernel install --user --name=viewer3d --display-name "Python (viewer3d)"

echo "✅ Setup complete!"
