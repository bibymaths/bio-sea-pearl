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

pretty_json() {
  jq .
}

write_fasta() {
  local path="$1"
  local header="$2"
  local seq="$3"
  printf ">%s\n%s\n" "$header" "$seq" > "$path"
}

read_file() {
  local path="$1"
  cat "$path"
}

post_json() {
  local endpoint="$1"
  local payload="$2"
  local body_file="$3"

  curl -sS \
    -o "$body_file" \
    -w "%{http_code}" \
    -X POST "${BASE_URL}${endpoint}" \
    -H "Content-Type: application/json" \
    -d "$payload"
}

get_json() {
  local endpoint="$1"
  local body_file="$2"

  curl -sS \
    -o "$body_file" \
    -w "%{http_code}" \
    "${BASE_URL}${endpoint}"
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
    echo "[PASS] $name response check OK: $expr"
  else
    echo "[FAIL] $name response check failed: $expr"
    echo "[DEBUG] Response body:"
    cat "$file" | pretty_json || cat "$file"
    return 1
  fi
}

show_payload_and_response() {
  local endpoint="$1"
  local payload="$2"
  local body_file="$3"

  echo "[INFO] POST $endpoint"
  echo "[INFO] Payload:"
  echo "$payload" | pretty_json || echo "$payload"
  echo "[INFO] Response body:"
  cat "$body_file" | pretty_json || cat "$body_file"
  echo
}

run_post_test() {
  local name="$1"
  local endpoint="$2"
  local payload="$3"
  local expected_status="$4"
  local jq_expr="$5"

  local body_file="$TMP_DIR/${name}.json"
  local status

  echo
  echo "===== $name ====="

  status="$(post_json "$endpoint" "$payload" "$body_file")"
  show_payload_and_response "$endpoint" "$payload" "$body_file"

  assert_status "$name" "$status" "$expected_status"
  assert_jq "$name" "$body_file" "$jq_expr"
}

run_get_test() {
  local name="$1"
  local endpoint="$2"
  local expected_status="$3"
  local jq_expr="$4"

  local body_file="$TMP_DIR/${name}.json"
  local status

  echo
  echo "===== $name ====="

  status="$(get_json "$endpoint" "$body_file")"

  echo "[INFO] GET $endpoint"
  echo "[INFO] Response body:"
  cat "$body_file" | pretty_json || cat "$body_file"
  echo

  assert_status "$name" "$status" "$expected_status"
  assert_jq "$name" "$body_file" "$jq_expr"
}

SEQ1_FA="$TMP_DIR/seq1.fa"
SEQ2_FA="$TMP_DIR/seq2.fa"
MARKOV_FA="$TMP_DIR/markov.fa"

write_fasta "$SEQ1_FA" "seq1" "ACGT"
write_fasta "$SEQ2_FA" "seq2" "ACCT"
write_fasta "$MARKOV_FA" "markov" "ACGTACGT"

SEQ1_CONTENT="$(read_file "$SEQ1_FA")"
SEQ2_CONTENT="$(read_file "$SEQ2_FA")"
MARKOV_CONTENT="$(read_file "$MARKOV_FA")"

run_get_test \
  "root" \
  "/" \
  "200" \
  '.message | type == "string"'

run_get_test \
  "openapi" \
  "/openapi.json" \
  "200" \
  '.paths["/align"] and .paths["/markov"] and .paths["/distance"] and .paths["/kmer"] and .paths["/bwt/search"]'

run_post_test \
  "distance" \
  "/distance" \
  '{"seq1":"kitten","seq2":"sitting","metric":"levenshtein"}' \
  "200" \
  '.distance == 3'

run_post_test \
  "kmer" \
  "/kmer" \
  '{"sequence":"ACGTACGT","k":3}' \
  "200" \
  '.counts["ACG"] == 2 and .counts["CGT"] == 2 and .counts["GTA"] == 1 and .counts["TAC"] == 1'

run_post_test \
  "bwt_search" \
  "/bwt/search" \
  '{"sequence":"ACGTACGT","pattern":"CGT"}' \
  "200" \
  '(.positions | type) == "array" and ((.positions | sort) == [1,5])'

ALIGN_PAYLOAD="$(jq -n \
  --arg fasta1 "$SEQ1_CONTENT" \
  --arg fasta2 "$SEQ2_CONTENT" \
  --arg matrix "alignment/scoring/blosum62.mat" \
  --arg mode "global" \
  '{fasta1:$fasta1,fasta2:$fasta2,matrix:$matrix,mode:$mode}')"

run_post_test \
  "align" \
  "/align" \
  "$ALIGN_PAYLOAD" \
  "200" \
  '.result | type == "string" and length > 0 and (contains("\u001b") | not)'

MARKOV_PAYLOAD="$(jq -n \
  --arg fasta "$MARKOV_CONTENT" \
  --argjson length 20 \
  --arg start "A" \
  --argjson order 1 \
  --arg method "alias" \
  --argjson pseudocount 0 \
  '{fasta:$fasta,length:$length,start:$start,order:$order,method:$method,pseudocount:$pseudocount}')"

run_post_test \
  "markov" \
  "/markov" \
  "$MARKOV_PAYLOAD" \
  "200" \
  '.walk | type == "string" and length > 0 and test("^[ACGT]+$")'

echo
echo "[DONE] All endpoint checks completed successfully."