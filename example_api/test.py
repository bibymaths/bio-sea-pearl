import requests
import json
import os
import sys

BASE_URL = "http://127.0.0.1:8000"


# ==========================================
# FILE READING HELPERS
# ==========================================
def read_fasta_file(filename):
    """Reads the entire file content as a string."""
    if not os.path.exists(filename):
        print(f"[ERROR] Could not find {filename} in the current directory.")
        sys.exit(1)
    with open(filename, "r") as file:
        return file.read()


def extract_sequence(fasta_content):
    """Strips the >header line and newlines to get the pure sequence string."""
    lines = fasta_content.strip().split('\n')
    # Join all lines that do NOT start with '>'
    sequence = "".join([line for line in lines if not line.startswith(">")])
    return sequence.strip()


def print_response(endpoint, response):
    print(f"\n--- Testing POST {endpoint} ---")
    if response.status_code == 200:
        print(json.dumps(response.json(), indent=2))
    else:
        print(f"Error {response.status_code}: {response.text}")


# ==========================================
# 1. ALIGNMENT (Using cleaned FASTA strings)
# ==========================================
def test_align(seq1_pure, seq2_pure):
    url = f"{BASE_URL}/align"

    # Create simple, safe FASTA formats to prevent legacy script crashes
    clean_fasta1 = f">seq1\n{seq1_pure}\n"
    clean_fasta2 = f">seq2\n{seq2_pure}\n"

    payload = {
        "fasta1": clean_fasta1,
        "fasta2": clean_fasta2,
        "matrix": "alignment/scoring/blosum62.mat",
        "mode": "global"
    }
    response = requests.post(url, json=payload)
    print_response("/align", response)


# ==========================================
# 2. MARKOV CHAIN (Using cleaned FASTA string)
# ==========================================
def test_markov(seq1_pure):
    url = f"{BASE_URL}/markov"

    clean_fasta = f">train\n{seq1_pure}\n"
    start_aa = seq1_pure[0] if seq1_pure else "M"

    payload = {
        "fasta": clean_fasta,
        "length": 20,
        "start": start_aa,
        "order": 1,
        "method": "alias",
        "pseudocount": 0
    }
    response = requests.post(url, json=payload)
    print_response("/markov", response)


# ==========================================
# 3. DISTANCE (Requires pure sequence strings)
# ==========================================
def test_distance(seq1_pure, seq2_pure):
    url = f"{BASE_URL}/distance"
    payload = {
        "seq1": seq1_pure,
        "seq2": seq2_pure,
        "metric": "levenshtein"
    }
    response = requests.post(url, json=payload)
    print_response("/distance", response)


# ==========================================
# 4. BWT / FM-INDEX SEARCH (Requires pure sequence string)
# ==========================================
def test_bwt_search(seq1_pure):
    url = f"{BASE_URL}/bwt/search"
    # We will search for the first 3 amino acids of seq1 so we are guaranteed a match
    search_pattern = seq1_pure[:3] if len(seq1_pure) >= 3 else seq1_pure

    payload = {
        "sequence": seq1_pure,
        "pattern": search_pattern
    }
    print(f"(Searching for pattern '{search_pattern}' in seq1)")
    response = requests.post(url, json=payload)
    print_response("/bwt/search", response)


# ==========================================
# 5. K-MER COUNTS (Requires pure sequence string)
# ==========================================
def test_kmer(seq1_pure):
    url = f"{BASE_URL}/kmer"
    payload = {
        "sequence": seq1_pure,
        "k": 2
    }
    response = requests.post(url, json=payload)
    print_response("/kmer", response)


# ==========================================
# MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    print(f"Connecting to Bio Sea Pearl API at {BASE_URL}...")

    # 1. Read the raw files
    print("Reading seq1.fa and seq2.fa from the root directory...")
    fasta1 = read_fasta_file("seq1.fa")
    fasta2 = read_fasta_file("seq2.fa")

    # 2. Extract the pure amino acid strings (remove headers)
    seq1_str = extract_sequence(fasta1)
    seq2_str = extract_sequence(fasta2)

    try:
        # 3. Run all tests passing the appropriate formats
        test_align(fasta1, fasta2)
        test_markov(seq1_str)
        test_distance(seq1_str, seq2_str)
        test_bwt_search(seq1_str)
        test_kmer(seq1_str)

    except requests.exceptions.ConnectionError:
        print(f"\n[ERROR] Could not connect to {BASE_URL}. Is your Uvicorn server running?")