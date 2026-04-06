#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "==> Building and starting Bio Sea Pearl containers..."
docker compose up --build -d

echo "==> Containers started.  API available at http://localhost:8000"
echo "    Logs: docker compose logs -f"
