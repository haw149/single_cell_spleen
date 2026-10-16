#!/bin/bash
set -euo pipefail
APP_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
if [ ! -x "$APP_DIR/.venv/bin/python" ]; then
  echo "Set up the environment first using the instructions in shiny_app/README.md."
  exit 1
fi
cd "$APP_DIR/.."
echo "EDNA will open at http://127.0.0.1:8765"
echo "Keep this terminal open. Press Control-C to stop the app."
exec "$APP_DIR/.venv/bin/python" -m shiny run --host 127.0.0.1 --port "${EDNA_PORT:-8765}" --launch-browser "$APP_DIR/app.py"
