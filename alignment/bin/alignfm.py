#!/usr/bin/env python3
import sys
import argparse
import os
from multiprocessing import Pool, Manager  # Manager won't be used in this revised approach
from functools import partial
import time  # For debugging or simple timeouts if needed

# Attempt to import FMIndex from the bwt module
try:
    # Assuming bwt.transform is the correct path if bwt is a package
    # If transform.py is in a directory named 'bwt' next to alignfm.py:
    from bwt.transform import FMIndex

    BWT_MODULE_AVAILABLE = True
except ImportError:
    # Fallback if bwt/transform.py is not found or if bwt is just transform.py
    try:
        from transform import FMIndex  # If transform.py is in the same dir or PYTHONPATH

        BWT_MODULE_AVAILABLE = True
    except ImportError:
        BWT_MODULE_AVAILABLE = False
        FMIndex = None  # Placeholder
        print("Warning: Could not import FMIndex. 'fm_local' mode will be unavailable.", file=sys.stderr)

# --- Global Configuration & Helper ---
GAP_OPEN = 10
GAP_EXTEND = 1
MODE = 'global'
MATRIX = {}
NEG_INF = -float('inf')

# For multiprocessing
IN_CHILD_PROCESS = False  # Initialize globally. True if in a Pool worker.
MAX_PROCS = int(os.environ.get('ALIGN_PROCS', 2))
# Heuristic for Hirschberg's NWScore, not directly for seed extensions here
MIN_SIZE_FOR_HIRSCHBERG_FORK = 50 * 50

# ANSI color constants
COL_GREEN = "\033[32m"
COL_CYAN = "\033[36m"
COL_RED = "\033[31m"
COL_RESET = "\033[0m"

_pool_globals = {}  # For worker initialization


def score_matrix_chars(char_a, char_b):
    # Ensure MATRIX is in _pool_globals, especially if called from main thread before pool init
    matrix_to_use = _pool_globals.get('MATRIX', MATRIX)  # Fallback to global MATRIX if not in _pool_globals
    return matrix_to_use.get(char_a, {}).get(char_b,
                                             matrix_to_use.get(char_b, {}).get(char_a, 0))


# --- File & Matrix I/O (largely unchanged, minor fixes) ---
def slurp_seq(filename):
    if filename == "-":
        fh = sys.stdin
    else:
        try:
            fh = open(filename, 'r')
        except FileNotFoundError:
            print(f"Error: Sequence file '{filename}' not found.", file=sys.stderr)
            sys.exit(1)

    seq_lines = []
    header_found = False
    first_seq_header = None
    for line in fh:
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            if first_seq_header is None:  # First header
                first_seq_header = line[1:]
                header_found = True
            else:  # Found another header, stop reading for this simple slurp
                break
            continue
        if header_found:  # Only append lines after the first header is found
            seq_lines.append(line)

    if filename != "-":
        fh.close()

    if not seq_lines:
        if filename == "-":
            print("Error: No sequence data read from stdin.", file=sys.stderr)
        else:
            print(f"Error: No sequence data read from '{filename}'. Is it a valid FASTA file with sequence data?",
                  file=sys.stderr)
        sys.exit(1)

    return "".join(seq_lines).upper()


def read_matrix_file(filename):
    mat = {}
    col_headers = []
    try:
        with open(filename, 'r') as fh:
            for line_num, line_content in enumerate(fh, 1):
                line = line_content.strip()
                if not line or line.startswith(('#', ';', '/', '*')):
                    continue

                tokens = [t for t in line.split() if t]
                if not tokens: continue

                if not col_headers:
                    if len(tokens) < 1:
                        print(f"Warning: Matrix header line '{line}' seems too short. Skipping.", file=sys.stderr)
                        continue
                    col_headers = tokens
                else:
                    if not tokens[0].isalpha() or len(tokens[0]) > 1:
                        print(
                            f"Warning: Skipping malformed matrix row {line_num}: '{line_content.strip()}'. Invalid row header.",
                            file=sys.stderr)
                        continue
                    if len(tokens) != len(col_headers) + 1:
                        print(
                            f"Warning: Skipping malformed matrix row {line_num}: '{line_content.strip()}'. Expected {len(col_headers) + 1} fields, got {len(tokens)}.",
                            file=sys.stderr)
                        continue

                    row_char = tokens[0]
                    mat[row_char] = {}
                    for i, col_char in enumerate(col_headers):
                        try:
                            mat[row_char][col_char] = int(tokens[i + 1])
                        except ValueError:
                            print(
                                f"Warning: Non-integer score for ({row_char}, {col_char}) in matrix row {line_num}. Using 0.",
                                file=sys.stderr)
                            mat[row_char][col_char] = 0
    except FileNotFoundError:
        print(f"Error: Scoring matrix file '{filename}' not found.", file=sys.stderr)
        sys.exit(1)

    if not col_headers:
        raise ValueError(f"No valid header line found in matrix '{filename}'")
    if not mat:
        raise ValueError(f"No score data found in matrix '{filename}'")
    return mat


# --- Core DP Algorithm (_NWScore_serial in Perl) ---
def gotoh_linear_space_scores(seq_a_list, seq_b_list, is_hirschberg_reversed_pass):
    # This function relies on _pool_globals being set correctly by init_worker or main thread
    m = len(seq_a_list)
    n = len(seq_b_list)

    current_seq_a = list(reversed(seq_a_list)) if is_hirschberg_reversed_pass else seq_a_list
    current_seq_b = list(reversed(seq_b_list)) if is_hirschberg_reversed_pass else seq_b_list

    s_prev = [0.0] * (n + 1)
    e_prev = [_pool_globals['NEG_INF']] * (n + 1)
    f_prev = [_pool_globals['NEG_INF']] * (n + 1)

    s_cur = [0.0] * (n + 1)
    e_cur = [0.0] * (n + 1)
    f_cur = [0.0] * (n + 1)

    s_prev[0] = 0.0
    e_prev[0], f_prev[0] = _pool_globals['NEG_INF'], _pool_globals['NEG_INF']

    for j in range(1, n + 1):
        if _pool_globals['MODE'] == 'global':
            s_prev[j] = -(_pool_globals['GAP_OPEN'] + (j - 1) * _pool_globals['GAP_EXTEND'])
            e_prev[j] = s_prev[j]
            f_prev[j] = _pool_globals['NEG_INF']
        elif _pool_globals['MODE'] == 'lcs':
            s_prev[j] = 0.0
            e_prev[j], f_prev[j] = _pool_globals['NEG_INF'], _pool_globals['NEG_INF']
        else:  # local
            s_prev[j], e_prev[j], f_prev[j] = 0.0, 0.0, 0.0

    for i in range(1, m + 1):
        if _pool_globals['MODE'] == 'global':
            s_cur[0] = -(_pool_globals['GAP_OPEN'] + (i - 1) * _pool_globals['GAP_EXTEND'])
            e_cur[0], f_cur[0] = _pool_globals['NEG_INF'], s_cur[0]
        elif _pool_globals['MODE'] == 'lcs':
            s_cur[0], e_cur[0], f_cur[0] = 0.0, _pool_globals['NEG_INF'], _pool_globals['NEG_INF']
        else:  # local
            s_cur[0], e_cur[0], f_cur[0] = 0.0, 0.0, 0.0

        for j in range(1, n + 1):
            char_a = current_seq_a[i - 1]
            char_b = current_seq_b[j - 1]

            sub_score = score_matrix_chars(char_a, char_b)
            if _pool_globals['MODE'] == 'lcs':
                sub_score = 1 if char_a == char_b else 0

            diag_pred_max = max(s_prev[j - 1], e_prev[j - 1], f_prev[j - 1])
            m_val = diag_pred_max + sub_score
            e_val = max(s_cur[j - 1] - _pool_globals['GAP_OPEN'], e_cur[j - 1] - _pool_globals['GAP_EXTEND'])
            f_val = max(s_prev[j] - _pool_globals['GAP_OPEN'], f_prev[j] - _pool_globals['GAP_EXTEND'])

            if _pool_globals['MODE'] == 'local':
                m_val, e_val, f_val = max(0, m_val), max(0, e_val), max(0, f_val)

            s_val = max(m_val, e_val, f_val)
            if _pool_globals['MODE'] == 'local': s_val = max(0, s_val)

            s_cur[j], e_cur[j], f_cur[j] = s_val, e_val, f_val

        s_prev, e_prev, f_prev = list(s_cur), list(e_cur), list(f_cur)

        # Progress bar logic: check global IN_CHILD_PROCESS
        if not IN_CHILD_PROCESS and not is_hirschberg_reversed_pass:
            progress_bar(i, m)

    return s_cur


# --- Multiprocessing Setup ---
_pool_instance = None  # Global pool instance


def init_worker(mode_init, gap_open_init, gap_extend_init, matrix_init, neg_inf_init):
    """Initializer for worker processes in the Pool."""
    global IN_CHILD_PROCESS  # Worker modifies the global IN_CHILD_PROCESS for its own context
    IN_CHILD_PROCESS = True
    _pool_globals['MODE'] = mode_init
    _pool_globals['GAP_OPEN'] = gap_open_init
    _pool_globals['GAP_EXTEND'] = gap_extend_init
    _pool_globals['MATRIX'] = matrix_init
    _pool_globals['NEG_INF'] = neg_inf_init


# --- Wrappers for multiprocessing Pool ---
def gotoh_linear_space_scores_wrapper(args_tuple):
    # This wrapper is called by a worker process. init_worker has set IN_CHILD_PROCESS=True.
    seq_a, seq_b, is_rev = args_tuple
    return gotoh_linear_space_scores(seq_a, seq_b, is_rev)


def small_align_core_wrapper(args_tuple):
    seq_a, seq_b = args_tuple
    return small_align_core(seq_a, seq_b)


def hirschberg_L_R_task_wrapper(args_tuple):
    seq_a_L, seq_b_L_common, rev_L, seq_a_R, rev_R = args_tuple
    l_scores = gotoh_linear_space_scores(seq_a_L, seq_b_L_common, rev_L)
    r_scores = gotoh_linear_space_scores(seq_a_R, seq_b_L_common, rev_R)
    return l_scores, r_scores


def extend_seed_hit_wrapper(args_tuple):
    """Wrapper for parallelizing seed extension."""
    query_sub_list, target_sub_list, orig_q_start, orig_t_start = args_tuple
    current_worker_mode = _pool_globals['MODE']
    try:
        _pool_globals['MODE'] = 'local'
        aln_q, aln_t = small_align_core(query_sub_list, target_sub_list)
    finally:
        _pool_globals['MODE'] = current_worker_mode

    score = 0
    gaps_q, gaps_t = 0, 0
    for i in range(len(aln_q)):
        cq, ct = aln_q[i], aln_t[i]
        if cq == '-':
            score -= (_pool_globals['GAP_OPEN'] if gaps_q == 0 else 0) + _pool_globals['GAP_EXTEND']
            gaps_q += 1
            gaps_t = 0
        elif ct == '-':
            score -= (_pool_globals['GAP_OPEN'] if gaps_t == 0 else 0) + _pool_globals['GAP_EXTEND']
            gaps_t += 1
            gaps_q = 0
        else:
            score += score_matrix_chars(cq, ct)
            gaps_q, gaps_t = 0, 0

    return aln_q, aln_t, score, orig_q_start, orig_t_start, len(query_sub_list), len(target_sub_list)


# --- Alignment Algorithms ---
def small_align_core(seq_a_list, seq_b_list):
    # This function relies on _pool_globals and global IN_CHILD_PROCESS
    m, n = len(seq_a_list), len(seq_b_list)
    s_mat = [[0.0] * (n + 1) for _ in range(m + 1)]
    e_mat = [[_pool_globals['NEG_INF']] * (n + 1) for _ in range(m + 1)]
    f_mat = [[_pool_globals['NEG_INF']] * (n + 1) for _ in range(m + 1)]
    p_mat = [[''] * (n + 1) for _ in range(m + 1)]

    s_mat[0][0], p_mat[0][0] = 0.0, 'X'

    for j in range(1, n + 1):
        if _pool_globals['MODE'] == 'global' or _pool_globals['MODE'] == 'lcs':
            val = -(_pool_globals['GAP_OPEN'] + (j - 1) * _pool_globals['GAP_EXTEND'])
            s_mat[0][j], e_mat[0][j], f_mat[0][j], p_mat[0][j] = val, val, _pool_globals['NEG_INF'], 'E'
        else:  # local
            s_mat[0][j], e_mat[0][j], f_mat[0][j], p_mat[0][j] = 0.0, 0.0, 0.0, 'X'

    for i in range(1, m + 1):
        if _pool_globals['MODE'] == 'global' or _pool_globals['MODE'] == 'lcs':
            val = -(_pool_globals['GAP_OPEN'] + (i - 1) * _pool_globals['GAP_EXTEND'])
            s_mat[i][0], e_mat[i][0], f_mat[i][0], p_mat[i][0] = val, _pool_globals['NEG_INF'], val, 'F'
        else:  # local
            s_mat[i][0], e_mat[i][0], f_mat[i][0], p_mat[i][0] = 0.0, 0.0, 0.0, 'X'

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            char_a, char_b = seq_a_list[i - 1], seq_b_list[j - 1]
            sub = score_matrix_chars(char_a, char_b)
            if _pool_globals['MODE'] == 'lcs': sub = 1 if char_a == char_b else 0

            diag_pred = max(s_mat[i - 1][j - 1], e_mat[i - 1][j - 1], f_mat[i - 1][j - 1])
            m_val_cand = diag_pred + sub
            e_val_cand = max(s_mat[i][j - 1] - _pool_globals['GAP_OPEN'], e_mat[i][j - 1] - _pool_globals['GAP_EXTEND'])
            f_val_cand = max(s_mat[i - 1][j] - _pool_globals['GAP_OPEN'], f_mat[i - 1][j] - _pool_globals['GAP_EXTEND'])

            if _pool_globals['MODE'] == 'local':
                m_val_cand, e_val_cand, f_val_cand = max(0, m_val_cand), max(0, e_val_cand), max(0, f_val_cand)

            e_mat[i][j], f_mat[i][j] = e_val_cand, f_val_cand

            best_score, ptr = m_val_cand, 'M'
            # Pointer preference: M > E > F (original Perl script might have implicit preferences)
            if e_val_cand > best_score:
                best_score, ptr = e_val_cand, 'E'
            # Check F even if E was better than M, to ensure F > E if that's the case
            if f_val_cand > best_score:
                best_score, ptr = f_val_cand, 'F'

            if _pool_globals['MODE'] == 'local' and 0 > best_score:
                s_mat[i][j], p_mat[i][j] = 0, 'X'
                e_mat[i][j], f_mat[i][j] = max(0, e_mat[i][j]), max(0, f_mat[i][j])
            else:
                s_mat[i][j], p_mat[i][j] = best_score, ptr

        if not IN_CHILD_PROCESS: progress_bar(i, m)  # Check global IN_CHILD_PROCESS

    align_a, align_b = "", ""
    cur_i, cur_j = m, n
    if _pool_globals['MODE'] == 'local':
        max_s = 0.0
        # Find the cell with the max score to start traceback for local alignment
        # Important: if multiple cells have max_score, this picks the one with largest i, then largest j.
        for r_idx_local in range(m + 1):
            for c_idx_local in range(n + 1):
                if s_mat[r_idx_local][c_idx_local] >= max_s:
                    max_s = s_mat[r_idx_local][c_idx_local]
                    cur_i, cur_j = r_idx_local, c_idx_local

    current_path = []  # For debugging traceback path
    while cur_i > 0 or cur_j > 0:
        current_path.append(p_mat[cur_i][cur_j])
        if p_mat[cur_i][cur_j] == 'X': break
        if p_mat[cur_i][cur_j] == 'M':
            align_a, align_b = seq_a_list[cur_i - 1] + align_a, seq_b_list[cur_j - 1] + align_b
            cur_i, cur_j = cur_i - 1, cur_j - 1
        elif p_mat[cur_i][cur_j] == 'E':  # Gap in A (came from left)
            align_a, align_b = '-' + align_a, seq_b_list[cur_j - 1] + align_b
            cur_j -= 1
        elif p_mat[cur_i][cur_j] == 'F':  # Gap in B (came from top)
            align_a, align_b = seq_a_list[cur_i - 1] + align_a, '-' + align_b
            cur_i -= 1
        else:  # Should be 'X' or error
            # print(f"Traceback anomaly: P[{cur_i}][{cur_j}] = {p_mat[cur_i][cur_j]} with path {current_path[::-1]}", file=sys.stderr)
            break
    return align_a, align_b


def hirschberg(seq_a_list, seq_b_list):
    # This function relies on _pool_globals and global IN_CHILD_PROCESS
    global IN_CHILD_PROCESS  # Allow modification for recursive calls' serial path
    m, n = len(seq_a_list), len(seq_b_list)

    if m == 0: return '-' * n, "".join(seq_b_list)
    if n == 0: return "".join(seq_a_list), '-' * m
    if m == 1 or n == 1:
        # Base case: use small_align_core.
        # Temporarily ensure serial execution if called from main thread.
        original_in_child_for_base = IN_CHILD_PROCESS
        try:
            IN_CHILD_PROCESS = True
            res = small_align_core(seq_a_list, seq_b_list)
        finally:
            IN_CHILD_PROCESS = original_in_child_for_base
        return res

    i_mid = m // 2
    seq_a_L_part = seq_a_list[0:i_mid]
    seq_a_R_part = seq_a_list[i_mid:m]

    args_for_L_R = (seq_a_L_part, seq_b_list, False, seq_a_R_part, True)

    should_fork_this_step = (not IN_CHILD_PROCESS and MAX_PROCS > 1 and
                             _pool_instance is not None and
                             m * n >= MIN_SIZE_FOR_HIRSCHBERG_FORK)

    l_score_row, r_score_row = None, None
    if should_fork_this_step:
        async_res = _pool_instance.apply_async(hirschberg_L_R_task_wrapper, args=(args_for_L_R,))
        try:
            res_tuple = async_res.get(timeout=7200)  # Consider making timeout configurable
            if res_tuple and len(res_tuple) == 2 and isinstance(res_tuple[0], list) and isinstance(res_tuple[1], list):
                l_score_row, r_score_row = res_tuple
            else:  # Fork returned unexpected data
                print(f"Hirschberg L/R fork for {m}x{n} returned invalid data. Falling back to serial.",
                      file=sys.stderr)
                l_score_row, r_score_row = None, None  # Ensure fallback is triggered
        except Exception as e:
            print(f"Hirschberg L/R fork error for {m}x{n}: {e}. Falling back to serial.", file=sys.stderr)
            l_score_row, r_score_row = None, None  # Ensure fallback is triggered

    if not l_score_row or not r_score_row:
        original_in_child_for_hirsch_serial = IN_CHILD_PROCESS
        try:
            IN_CHILD_PROCESS = True  # Ensure serial execution for these gotoh calls
            l_score_row = gotoh_linear_space_scores(seq_a_L_part, seq_b_list, False)
            r_score_row = gotoh_linear_space_scores(seq_a_R_part, seq_b_list, True)
        finally:
            IN_CHILD_PROCESS = original_in_child_for_hirsch_serial

    if not l_score_row or not r_score_row:  # Should not happen if gotoh works
        raise Exception(
            f"Critical error: Failed to compute score arrays for Hirschberg {m}x{n} even in serial fallback.")

    best_score_sum, j_mid = _pool_globals['NEG_INF'], 0
    for j in range(n + 1):  # Iterate 0 to n (inclusive for B_j)
        score_l = l_score_row[j] if j < len(l_score_row) else _pool_globals['NEG_INF']
        # r_score_row is for reversed B, so index is (n-j)
        score_r_idx = n - j
        score_r = r_score_row[score_r_idx] if score_r_idx < len(r_score_row) else _pool_globals['NEG_INF']

        current_sum = score_l + score_r
        if current_sum > best_score_sum:
            best_score_sum, j_mid = current_sum, j

    # Recursive calls will manage their own IN_CHILD_PROCESS state if they fork
    a1, b1 = hirschberg(seq_a_list[0:i_mid], seq_b_list[0:j_mid])
    a2, b2 = hirschberg(seq_a_list[i_mid:m], seq_b_list[j_mid:n])  # j_mid can be n
    return a1 + a2, b1 + b2


def seed_and_extend_alignment(query_seq_list, target_seq_list, fm_idx_target,
                              seed_len=12, seed_step=1,
                              extension_window_size=100):
    """
    Performs local alignment using FM-index for seeding and Smith-Waterman for extension.
    """
    global IN_CHILD_PROCESS  # This function reads the global IN_CHILD_PROCESS

    if not BWT_MODULE_AVAILABLE:
        print("Error: BWT module (for FMIndex) is not available. Cannot use 'fm_local' mode.", file=sys.stderr)
        return "", "", _pool_globals['NEG_INF']

    query_str = "".join(query_seq_list)
    best_overall_score = _pool_globals['NEG_INF']  # Or 0.0 for local if no positive score found
    best_align_q, best_align_t = "", ""
    seed_extension_tasks = []

    print(f"Seeding: query_len={len(query_seq_list)}, target_len={len(target_seq_list)}, seed_len={seed_len}",
          file=sys.stderr)
    num_seeds_processed = 0
    for i in range(0, len(query_str) - seed_len + 1, seed_step):
        seed = query_str[i: i + seed_len]
        num_seeds_processed += 1

        l, r = fm_idx_target.backward_search(seed)
        if l < r:
            target_hits_sa_indices = fm_idx_target.sa[l:r]
            for hit_sa_idx in target_hits_sa_indices:
                if hit_sa_idx >= len(target_seq_list): continue

                q_start_seed = i
                t_start_seed = hit_sa_idx

                q_win_start = max(0, q_start_seed - extension_window_size // 2)
                q_win_end = min(len(query_seq_list), q_start_seed + seed_len + extension_window_size // 2)
                sub_query_list = query_seq_list[q_win_start: q_win_end]

                t_win_start = max(0, t_start_seed - extension_window_size // 2)
                t_win_end = min(len(target_seq_list), t_start_seed + seed_len + extension_window_size // 2)
                sub_target_list = target_seq_list[t_win_start: t_win_end]

                if not sub_query_list or not sub_target_list: continue
                task_args = (sub_query_list, sub_target_list, q_win_start, t_win_start)
                seed_extension_tasks.append(task_args)

    print(f"Generated {len(seed_extension_tasks)} seed extension tasks from {num_seeds_processed} seeds.",
          file=sys.stderr)

    if not seed_extension_tasks:
        print("No seed hits found or tasks generated.", file=sys.stderr)
        return "", "", 0.0

    results = []
    # Check global IN_CHILD_PROCESS here
    if not IN_CHILD_PROCESS and MAX_PROCS > 1 and _pool_instance is not None and seed_extension_tasks:
        print(f"Extending {len(seed_extension_tasks)} seeds in parallel ({MAX_PROCS} processes)...", file=sys.stderr)
        map_results = _pool_instance.map_async(extend_seed_hit_wrapper, seed_extension_tasks)

        # Progress for map_async
        total_tasks = len(seed_extension_tasks)
        # Note: map_async doesn't have a direct way to get # completed tasks easily.
        # This loop is just a wait with a timeout.
        # For real progress, one might use `apply_async` with callbacks or `imap_unordered` and count.
        start_time = time.time()
        while not map_results.ready():
            time.sleep(0.5)  # Wait a bit
            if time.time() - start_time > 7200:  # 2hr timeout for all extensions
                print("\nTimeout waiting for seed extensions.", file=sys.stderr)
                _pool_instance.terminate()  # Forcefully stop pool on timeout
                _pool_instance.join()
                return best_align_q, best_align_t, best_overall_score  # Return what we have

        print("\nAll seed extensions completed (parallel).", file=sys.stderr)
        try:
            results = map_results.get()
        except Exception as e:
            print(f"Error getting results from parallel seed extensions: {e}", file=sys.stderr)
            # Potentially partial results might be available if some tasks succeeded before error
            # but map_results.get() will raise if any task raised an unhandled exception.

    else:
        print(f"Extending {len(seed_extension_tasks)} seeds serially...", file=sys.stderr)
        original_in_child_for_serial_extend = IN_CHILD_PROCESS
        try:
            IN_CHILD_PROCESS = True  # Ensure small_align_core runs serially within this loop
            for idx, task_args in enumerate(seed_extension_tasks):
                results.append(extend_seed_hit_wrapper(task_args))
                progress_bar(idx + 1, len(seed_extension_tasks))
        finally:
            IN_CHILD_PROCESS = original_in_child_for_serial_extend

    current_best_local_score = 0.0  # Local alignments start from 0
    final_aln_q, final_aln_t = "", ""

    for result in results:
        if result is None: continue
        aln_q_sub, aln_t_sub, score_sub, _, _, _, _ = result
        if score_sub > current_best_local_score:
            current_best_local_score = score_sub
            final_aln_q, final_aln_t = aln_q_sub, aln_t_sub

    return final_aln_q, final_aln_t, current_best_local_score


# --- Output and Utility Functions ---
def extract_lcs_from_alignment(aln_a, aln_b):
    runs, current_run = [], ""
    for char_a, char_b in zip(aln_a, aln_b):
        if char_a == char_b and char_a != '-':
            current_run += char_a
        else:
            if current_run: runs.append(current_run)
            current_run = ""
    if current_run: runs.append(current_run)
    return runs


def build_full_dp_matrix(seq_a_list, seq_b_list):
    global MODE  # Ensure MODE is accessible, should be from _pool_globals or main's args
    # This function uses _pool_globals or script globals if not in worker
    mode_to_use = _pool_globals.get('MODE', MODE)
    gap_open_to_use = _pool_globals.get('GAP_OPEN', GAP_OPEN)
    gap_extend_to_use = _pool_globals.get('GAP_EXTEND', GAP_EXTEND)
    neg_inf_to_use = _pool_globals.get('NEG_INF', NEG_INF)

    m, n = len(seq_a_list), len(seq_b_list)
    full_dp_S = [[0.0] * (n + 1) for _ in range(m + 1)]
    s_prev = [0.0] * (n + 1)
    e_prev = [neg_inf_to_use] * (n + 1)
    f_prev = [neg_inf_to_use] * (n + 1)

    s_cur = [0.0] * (n + 1)  # Not strictly needed here as we only fill full_dp_S
    e_cur = [0.0] * (n + 1)
    f_cur = [0.0] * (n + 1)

    s_prev[0], e_prev[0], f_prev[0] = 0.0, neg_inf_to_use, neg_inf_to_use
    for j in range(1, n + 1):
        if mode_to_use == 'global':
            s_prev[j] = -(gap_open_to_use + (j - 1) * gap_extend_to_use)
            e_prev[j], f_prev[j] = s_prev[j], neg_inf_to_use
        elif mode_to_use == 'lcs':
            s_prev[j], e_prev[j], f_prev[j] = 0.0, neg_inf_to_use, neg_inf_to_use
        else:  # local
            s_prev[j], e_prev[j], f_prev[j] = 0.0, 0.0, 0.0
    full_dp_S[0] = list(s_prev)

    for i in range(1, m + 1):
        # Initialize column 0 for current row i
        s_val_col0 = 0.0
        e_val_col0 = 0.0
        f_val_col0 = 0.0

        if mode_to_use == 'global':
            s_val_col0 = -(gap_open_to_use + (i - 1) * gap_extend_to_use)
            e_val_col0 = neg_inf_to_use
            f_val_col0 = s_val_col0
        elif mode_to_use == 'lcs':
            s_val_col0 = 0.0
            e_val_col0, f_val_col0 = neg_inf_to_use, neg_inf_to_use
        else:  # local
            s_val_col0, e_val_col0, f_val_col0 = 0.0, 0.0, 0.0

        # Set S[i][0] in full_dp_S
        full_dp_S[i][0] = s_val_col0
        # Update e_cur[0] and f_cur[0] for the DP calculation of the first cell S[i][1]
        e_cur[0] = e_val_col0
        f_cur[0] = f_val_col0

        for j in range(1, n + 1):
            char_a, char_b = seq_a_list[i - 1], seq_b_list[j - 1]
            sub = score_matrix_chars(char_a, char_b)
            if mode_to_use == 'lcs': sub = 1 if char_a == char_b else 0

            # S_prev refers to full_dp_S[i-1]
            # S_cur for [i][j-1] refers to full_dp_S[i][j-1]
            diag_pred = max(full_dp_S[i - 1][j - 1], e_prev[j - 1], f_prev[j - 1])
            m_val = diag_pred + sub

            e_val = max(full_dp_S[i][j - 1] - gap_open_to_use, e_cur[j - 1] - gap_extend_to_use)
            f_val = max(full_dp_S[i - 1][j] - gap_open_to_use, f_prev[j] - gap_extend_to_use)

            if mode_to_use == 'local':
                m_val, e_val, f_val = max(0, m_val), max(0, e_val), max(0, f_val)

            s_val = max(m_val, e_val, f_val)
            if mode_to_use == 'local': s_val = max(0, s_val)

            full_dp_S[i][j] = s_val
            e_cur[j] = e_val  # Store for next iteration's e_prev
            f_cur[j] = f_val  # Store for next iteration's f_prev

        # After row i is computed, its e_cur and f_cur become e_prev and f_prev for row i+1
        e_prev = list(e_cur)
        f_prev = list(f_cur)
        # s_prev for the next row (i+1) will be full_dp_S[i]

    return full_dp_S


def print_colored_alignment(aln_a, aln_b):
    len_aln = len(aln_a)
    if len_aln == 0: sys.stdout.write("No alignment to print.\n"); return
    pos_a, pos_b = 0, 0
    for offset in range(0, len_aln, 80):
        block_a, block_b = aln_a[offset: offset + 80], aln_b[offset: offset + 80]
        blen = len(block_a)
        sys.stdout.write(f"{pos_a + 1:6d} ")
        current_line_pos_a = pos_a
        for k in range(blen):
            ca, cb = block_a[k], block_b[k]
            color = COL_RED if ca == '-' or cb == '-' else COL_GREEN if ca == cb else COL_CYAN
            sys.stdout.write(f"{color}{ca}{COL_RESET}")
            if ca != '-': current_line_pos_a += 1
        sys.stdout.write(f" {current_line_pos_a}\n")
        sys.stdout.write(f"{pos_b + 1:6d} ")
        current_line_pos_b = pos_b
        for k in range(blen):
            ca, cb = block_a[k], block_b[k]
            color = COL_RED if ca == '-' or cb == '-' else COL_GREEN if ca == cb else COL_CYAN
            sys.stdout.write(f"{color}{cb}{COL_RESET}")
            if cb != '-': current_line_pos_b += 1
        sys.stdout.write(f" {current_line_pos_b}\n\n")
        pos_a, pos_b = current_line_pos_a, current_line_pos_b
    sys.stdout.flush()


def progress_bar(done, total, width=40):
    if total == 0: return
    pct = min(1, max(0, done / total))
    filled = int(pct * width)
    bar = "[" + "#" * filled + " " * (width - filled) + "]"
    sys.stderr.write(f"\r{bar} {pct * 100:3.0f}%")
    if done >= total: sys.stderr.write("\n")
    sys.stderr.flush()


# --- Main Application Logic ---
def main():
    global MAX_PROCS, _pool_instance, _pool_globals, NEG_INF, IN_CHILD_PROCESS, MODE, GAP_OPEN, GAP_EXTEND, MATRIX

    _pool_globals['NEG_INF'] = NEG_INF

    parser = argparse.ArgumentParser(description="Pairwise sequence alignment tool (Python port).")
    parser.add_argument("seq_a_file", help="Input FASTA file for sequence A")
    parser.add_argument("seq_b_file", help="Input FASTA file for sequence B")
    parser.add_argument("--mode", choices=['global', 'local', 'lcs', 'fm_local'], default='global',
                        help="Alignment mode (default: global; fm_local uses seed-and-extend)")
    parser.add_argument("--matrix", required=True, help="Scoring matrix file")
    parser.add_argument("--gapopen", type=int, default=10, help="Gap open penalty (positive value)")
    parser.add_argument("--gapext", type=int, default=1, help="Gap extend penalty (positive value)")
    parser.add_argument("--out", default="align", help="Output prefix for files")
    parser.add_argument("--maxprocs", type=int, default=MAX_PROCS,
                        help=f"Max processes for parallel parts (default: {MAX_PROCS})")
    parser.add_argument("--sentinel", default="$", help="Sentinel character for FM-index (if using fm_local)")
    parser.add_argument("--seedlen", type=int, default=2, help="Seed length for fm_local mode")
    parser.add_argument("--seedstep", type=int, default=1, help="Step for selecting seeds in fm_local mode")
    parser.add_argument("--extwindow", type=int, default=100, help="Extension window size around seeds for fm_local")

    args = parser.parse_args()

    # Set script-level globals from args, and also _pool_globals for workers and direct calls
    MODE = args.mode
    GAP_OPEN = args.gapopen
    GAP_EXTEND = args.gapext
    MAX_PROCS = args.maxprocs  # Update global MAX_PROCS

    _pool_globals['MODE'] = MODE
    _pool_globals['GAP_OPEN'] = GAP_OPEN
    _pool_globals['GAP_EXTEND'] = GAP_EXTEND

    try:
        MATRIX = read_matrix_file(args.matrix)
        _pool_globals['MATRIX'] = MATRIX
    except Exception as e:
        print(f"Error reading matrix file '{args.matrix}': {e}", file=sys.stderr)
        sys.exit(1)

    seq_a_str = slurp_seq(args.seq_a_file)
    seq_b_str = slurp_seq(args.seq_b_file)
    if not seq_a_str or not seq_b_str:
        print("Error: One or both input sequences are empty.", file=sys.stderr)
        sys.exit(1)

    seq_a_list = list(seq_a_str)
    seq_b_list = list(seq_b_str)

    if MAX_PROCS > 1:
        init_args_tuple = (
            _pool_globals['MODE'], _pool_globals['GAP_OPEN'], _pool_globals['GAP_EXTEND'],
            _pool_globals['MATRIX'], _pool_globals['NEG_INF']
        )
        # Important: set maxtasksperchild to prevent memory leaks in long-running workers
        _pool_instance = Pool(processes=MAX_PROCS, initializer=init_worker, initargs=init_args_tuple,
                              maxtasksperchild=100)
        # IN_CHILD_PROCESS is False by default for the main process

    aligned_s1, aligned_s2 = "", ""
    final_score_val = _pool_globals['NEG_INF']

    if args.mode == 'fm_local':
        if not BWT_MODULE_AVAILABLE:
            print("Error: BWT module required for 'fm_local' mode is not installed or found.", file=sys.stderr)
            if _pool_instance: _pool_instance.close(); _pool_instance.join()
            sys.exit(1)
        if args.sentinel in seq_b_str or args.sentinel in seq_a_str:
            print(f"Error: Sentinel character '{args.sentinel}' found in input sequences.", file=sys.stderr)
            if _pool_instance: _pool_instance.close(); _pool_instance.join()
            sys.exit(1)

        print(f"Building FM-index for sequence B (target, length {len(seq_b_str)})...", file=sys.stderr)
        fm_idx_b = FMIndex(seq_b_str, sentinel=args.sentinel)
        print("FM-index built.", file=sys.stderr)

        aligned_s1, aligned_s2, final_score_val = seed_and_extend_alignment(
            seq_a_list, seq_b_list, fm_idx_b,
            seed_len=args.seedlen, seed_step=args.seedstep,
            extension_window_size=args.extwindow
        )

    elif args.mode == 'lcs':  # Use args.mode for condition, _pool_globals['MODE'] for algorithm logic
        original_in_child_for_lcs = IN_CHILD_PROCESS
        try:
            IN_CHILD_PROCESS = True
            aligned_s1, aligned_s2 = small_align_core(seq_a_list, seq_b_list)
        finally:
            IN_CHILD_PROCESS = original_in_child_for_lcs
    else:  # global or local (non-FM)
        aligned_s1, aligned_s2 = hirschberg(seq_a_list, seq_b_list)

    if _pool_instance:
        _pool_instance.close()
        _pool_instance.join()

    print("\n")
    print_colored_alignment(aligned_s1, aligned_s2)

    with open(f"{args.out}.A.fa", "w") as f:
        f.write(f">A_aligned_{args.mode}\n{aligned_s1}\n")
    with open(f"{args.out}.B.fa", "w") as f:
        f.write(f">B_aligned_{args.mode}\n{aligned_s2}\n")
    print(f"Aligned sequences saved to {args.out}.A.fa and {args.out}.B.fa", file=sys.stderr)

    if args.mode != 'fm_local':
        # For non-fm_local modes, calculate score using full DP matrix
        # Ensure build_full_dp_matrix uses the correct mode for scoring
        # It will use the mode set in _pool_globals, which matches args.mode here.
        dp_full = build_full_dp_matrix(seq_a_list, seq_b_list)
        if dp_full:
            m_dp, n_dp = len(dp_full) - 1, (len(dp_full[0]) - 1 if dp_full and dp_full[0] else -1)
            if n_dp >= 0:  # Check if matrix has columns
                if args.mode == 'local':  # Use args.mode for this logic
                    current_max_score = 0.0  # Local score cannot be less than 0
                    for r_idx in range(m_dp + 1):
                        for c_idx in range(n_dp + 1):
                            if dp_full[r_idx][c_idx] > current_max_score:
                                current_max_score = dp_full[r_idx][c_idx]
                    final_score_val = current_max_score
                else:  # global or lcs
                    final_score_val = dp_full[m_dp][n_dp]
            else:
                final_score_val = "Error: DP matrix malformed"
        else:
            final_score_val = "Error: DP matrix not built"

        if dp_full:
            try:
                with open(f"{args.out}.matrix.tsv", "w") as f_tsv:
                    f_tsv.write("\t" + "\t".join(seq_b_list) + "\n")
                    for i, row_scores in enumerate(dp_full):
                        label = seq_a_list[i - 1] if i > 0 else ''
                        f_tsv.write(label + "\t" + "\t".join(map(str, row_scores)) + "\n")
                print(f"Full DP matrix saved to {args.out}.matrix.tsv", file=sys.stderr)
            except Exception as e:
                print(f"Error writing TSV matrix: {e}", file=sys.stderr)

    if args.mode == 'lcs':
        lcs_runs = extract_lcs_from_alignment(aligned_s1, aligned_s2)
        lcs_out_file = f"{args.out}.lcs.fa"
        if lcs_runs:
            longest_lcs = max(lcs_runs, key=len) if lcs_runs else ""
            with open(lcs_out_file, "w") as f_lcs:
                f_lcs.write(f">{args.out}_LCS\n{longest_lcs}\n")
            print(f"Longest LCS saved to {lcs_out_file}", file=sys.stderr)
        else:
            with open(lcs_out_file, "w") as f_lcs:
                f_lcs.write(f">{args.out}_LCS\n\n")
            print(f"No common subsequence found for {lcs_out_file}", file=sys.stderr)

    print(f"\nScore   :  {final_score_val}")
    print(f"Output  :  {args.out}.A.fa, {args.out}.B.fa")
    print(f"Completed {args.mode} alignment.")


if __name__ == "__main__":
    main()