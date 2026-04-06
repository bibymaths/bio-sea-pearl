#!/usr/bin/env bash
set -euo pipefail

BASE_URL="${BASE_URL:-http://127.0.0.1:8000}"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

echo "[INFO] Base URL: $BASE_URL"
echo "[INFO] Temp dir: $TMP_DIR"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "[ERROR] Required command not found: $1" >&2
    exit 1
  }
}

require_cmd curl
require_cmd jq

# Optional helpers
HAS_PYTHON=0
if command -v python >/dev/null 2>&1; then
  HAS_PYTHON=1
fi

write_fasta() {
  local path="$1"
  local header="$2"
  local seq="$3"
  printf ">%s\n%s\n" "$header" "$seq" > "$path"
}

SEQ1_FA="$TMP_DIR/seq1.fa"
SEQ2_FA="$TMP_DIR/seq2.fa"
MARKOV_FA="$TMP_DIR/markov.fa"

write_fasta "$SEQ1_FA" "seq1" "ACGT"
write_fasta "$SEQ2_FA" "seq2" "ACCT"
write_fasta "$MARKOV_FA" "markov" "ACGTACGT"

pretty_json() {
  jq .
}

post_json() {
  local endpoint="$1"
  local payload="$2"
  local body_file="$3"
  local status

  status="$(
    curl -sS \
      -o "$body_file" \
      -w "%{http_code}" \
      -X POST "${BASE_URL}${endpoint}" \
      -H "Content-Type: application/json" \
      -d "$payload"
  )"

  echo "$status"
}

assert_status() {
  local name="$1"
  local actual="$2"
  local expected="$3"

  if [[ "$actual" != "$expected" ]]; then
    echo "[FAIL] $name returned HTTP $actual, expected $expected"
    return 1
  fi
  echo "[PASS] $name returned HTTP $actual"
}

assert_jq() {
  local name="$1"
  local file="$2"
  local expr="$3"

  if jq -e "$expr" "$file" >/dev/null; then
    echo "[PASS] $name response shape OK: $expr"
  else
    echo "[FAIL] $name response shape check failed: $expr"
    echo "[DEBUG] Response body:"
    cat "$file" | pretty_json || cat "$file"
    return 1
  fi
}

run_test() {
  local name="$1"
  local endpoint="$2"
  local payload="$3"
  local expected_status="$4"
  local jq_expr="$5"

  local body_file="$TMP_DIR/${name}.json"

  echo
  echo "===== $name ====="
  echo "[INFO] POST $endpoint"
  echo "[INFO] Payload:"
  echo "$payload" | pretty_json || echo "$payload"

  local status
  status="$(post_json "$endpoint" "$payload" "$body_file")"

  echo "[INFO] Response body:"
  cat "$body_file" | pretty_json || cat "$body_file"
  echo

  assert_status "$name" "$status" "$expected_status"
  assert_jq "$name" "$body_file" "$jq_expr"
}

check_root() {
  echo
  echo "===== root ====="
  local body_file="$TMP_DIR/root.json"
  local status

  status="$(
    curl -sS \
      -o "$body_file" \
      -w "%{http_code}" \
      "${BASE_URL}/"
  )"

  echo "[INFO] GET /"
  echo "[INFO] Response body:"
  cat "$body_file" | pretty_json || cat "$body_file"
  echo

  assert_status "root" "$status" "200"

  if jq -e 'type == "object" or type == "string"' "$body_file" >/dev/null 2>&1; then
    echo "[PASS] root response is valid JSON"
  else
    echo "[FAIL] root response is not valid JSON"
    return 1
  fi
}

check_openapi() {
  echo
  echo "===== openapi ====="
  local body_file="$TMP_DIR/openapi.json"
  local status

  status="$(
    curl -sS \
      -o "$body_file" \
      -w "%{http_code}" \
      "${BASE_URL}/openapi.json"
  )"

  assert_status "openapi" "$status" "200"
  assert_jq "openapi" "$body_file" '
    .paths["/align"] and
    .paths["/markov"] and
    .paths["/distance"] and
    .paths["/kmer"] and
    .paths["/bwt/search"]
  '
}

check_root
check_openapi

run_test \
  "distance" \
  "/distance" \
  '{"seq1":"kitten","seq2":"sitting","metric":"levenshtein"}' \
  "200" \
  '.distance == 3'

run_test \
  "kmer" \
  "/kmer" \
  '{"sequence":"ACGTACGT","k":3}' \
  "200" \
  '.counts["ACG"] == 2 and .counts["CGT"] == 2 and .counts["GTA"] == 1 and .counts["TAC"] == 1'

run_test \
  "bwt_search" \
  "/bwt/search" \
  '{"sequence":"ACGTACGT","pattern":"CGT"}' \
  "200" \
  '(.positions | type) == "array" and (.positions | length) == 2'

ALIGN_PAYLOAD="$(jq -n \
  --arg fasta1 "$SEQ1_FA" \
  --arg fasta2 "$SEQ2_FA" \
  --arg matrix "alignment/scoring/blosum62.mat" \
  --arg mode "global" \
  '{fasta1:$fasta1,fasta2:$fasta2,matrix:$matrix,mode:$mode}')"

run_test \
  "align" \
  "/align" \
  "$ALIGN_PAYLOAD" \
  "200" \
  'type == "object"'

MARKOV_PAYLOAD="$(jq -n \
  --arg fasta "$MARKOV_FA" \
  --argjson length 20 \
  --arg start "A" \
  --argjson order 1 \
  --arg method "alias" \
  --argjson pseudocount 0 \
  '{fasta:$fasta,length:$length,start:$start,order:$order,method:$method,pseudocount:$pseudocount}')"

run_test \
  "markov" \
  "/markov" \
  "$MARKOV_PAYLOAD" \
  "200" \
  'type == "object"'

echo
echo "[DONE] All endpoint checks completed."