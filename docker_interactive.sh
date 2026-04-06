#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CONTAINER_NAME="biosea"

# If the container is already running, exec into it; otherwise start one.
if docker compose ps --status running | grep -q "$CONTAINER_NAME"; then
    echo "==> Opening shell in running container..."
    docker compose exec "$CONTAINER_NAME" /bin/bash
else
    echo "==> No running container found. Starting a temporary one..."
    docker compose run --rm "$CONTAINER_NAME" /bin/bash
fi
