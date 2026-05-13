#!/usr/bin/env bash
set -e

ROOT="$(cd "$(dirname "$0")" && pwd)"

# S'assure que npm est dans le PATH (installation Homebrew)
export PATH="/opt/homebrew/bin:$PATH"

cleanup() {
  echo ""
  echo "Arrêt..."
  kill "$API_PID" "$FRONT_PID" 2>/dev/null
  wait "$API_PID" "$FRONT_PID" 2>/dev/null
  exit 0
}
trap cleanup INT TERM

echo "▶ API      → http://localhost:8000"
cd "$ROOT"
uv run uvicorn api.main:app --host 0.0.0.0 --port 8000 &
API_PID=$!

echo "▶ Frontend → http://localhost:5173"
cd "$ROOT/frontend"
npm run dev &
FRONT_PID=$!

echo ""
echo "Ctrl+C pour tout arrêter"
wait
