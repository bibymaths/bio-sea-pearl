#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "==> Stopping and removing Bio Sea Pearl containers..."
docker compose down --remove-orphans

echo "==> Containers stopped."
