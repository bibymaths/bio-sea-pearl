#!/usr/bin/env python3
import sys
import argparse
import os
from multiprocessing import Pool, Manager
from functools import partial

# --- Global Configuration & Helper ---
# Will be set in main based on args
GAP_OPEN = 10
GAP_EXTEND = 1
MODE = 'global'  # global, local, lcs
MATRIX = {}
NEG_INF = -float('inf')

# For multiprocessing
IN_CHILD = False  # Emulate Perl's $IN_CHILD behavior
MAX_PROCS = 16
MIN_SIZE_FOR_FORK = 50 * 50  # Heuristic: m*n, don't fork for smaller problems

# ANSI color constants
COL_GREEN = "\033[32m"
COL_CYAN = "\033[36m"
COL_RED = "\033[31m"
COL_RESET = "\033[0m"

# Store forked results globally, managed by multiprocessing.Manager
# This will be initialized in main if MAX_PROCS > 0
_fork_results_manager = None
FORK_RESULTS = None


def score_matrix_chars(char_a, char_b):
    """Gets score from the global MATRIX."""
    return MATRIX.get(char_a, {}).get(char_b, MATRIX.get(char_b, {}).get(char_a, 0))


# --- File & Matrix I/O ---
def slurp_seq(filename):
    """Reads a sequence from a FASTA file, removing headers and whitespace."""
    if filename == "-":
        fh = sys.stdin
    else:
        fh = open(filename, 'r')

    seq_lines = []
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if not line.startswith('>'):
            seq_lines.append(line)

    if filename != "-":
        fh.close()
    return "".join(seq_lines).upper()


def read_matrix_file(filename):
    """Reads a scoring matrix from a file."""
    mat = {}
    col_headers = []
    with open(filename, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(('#', ';', '/', '*')):
                continue

            tokens = [t for t in line.split() if t]
            if not tokens:
                continue

            if not col_headers:  # First non-comment line is column headers
                if len(tokens) < 2:
                    raise ValueError("Matrix header line must have at least two characters.")
                col_headers = tokens
            else:  # Score lines
                if len(tokens) != len(col_headers) + 1:
                    # Allow for matrices that only define the upper/lower triangle
                    # For simplicity now, expect full definition as in Perl script
                    # raise ValueError(f"Matrix row '{line}' has incorrect number of fields. Expected {len(col_headers)+1}, got {len(tokens)}")
                    pass  # Silently skip malformed rows for now, or handle more gracefully

                row_char = tokens[0]
                mat[row_char] = {}
                for i, col_char in enumerate(col_headers):
                    if i + 1 < len(tokens):  # Check if token exists
                        try:
                            mat[row_char][col_char] = int(tokens[i + 1])
                        except ValueError:
                            # Handle cases where a score might not be an int, e.g. '*'
                            # For now, assume int scores or use a default.
                            mat[row_char][col_char] = 0  # Default for non-integer scores

    if not col_headers:
        raise ValueError(f"No valid header line found in matrix '{filename}'")
    if not mat:
        raise ValueError(f"No score data found in matrix '{filename}'")
    return mat


# --- Core DP Algorithm (_NWScore_serial in Perl) ---
def gotoh_linear_space_scores(seq_a_list, seq_b_list, is_reversed_for_hirschberg):
    """
    Computes scores for Gotoh algorithm in linear space.
    This corresponds to _NWScore_serial in the Perl script.
    Returns the last row of scores (S_cur).
    """
    global IN_CHILD  # To control progress bar display

    m = len(seq_a_list)
    n = len(seq_b_list)

    if is_reversed_for_hirschberg:
        # For Hirschberg, the actual sequence reversal is handled by slicing.
        # This flag is more about which half it is.
        pass

    # S, E, F arrays for previous and current rows
    s_prev = [0.0] * (n + 1)
    e_prev = [NEG_INF] * (n + 1)  # e_prev is not really used for row 0 init of S
    f_prev = [NEG_INF] * (n + 1)  # f_prev is not really used for row 0 init of S

    s_cur = [0.0] * (n + 1)
    e_cur = [0.0] * (n + 1)  # Stores E_ij (gap in seq_a_list)
    f_cur = [0.0] * (n + 1)  # Stores F_ij (gap in seq_b_list)

    # Initialize row 0 (i=0)
    s_prev[0] = 0.0
    e_prev[0] = NEG_INF  # E[0,0]
    f_prev[0] = NEG_INF  # F[0,0]

    for j in range(1, n + 1):
        if MODE == 'global':
            # Cost of j gaps in seq_a_list when aligning to empty prefix of seq_b_list
            s_prev[j] = -(GAP_OPEN + (j - 1) * GAP_EXTEND)
            e_prev[j] = s_prev[j]  # E[0,j] ends with gap in A, filled from S[0,j-1] - GO or E[0,j-1]-GE
            f_prev[j] = NEG_INF
        elif MODE == 'lcs':
            s_prev[j] = 0.0
            e_prev[j] = NEG_INF
            f_prev[j] = NEG_INF
        else:  # local
            s_prev[j] = 0.0
            e_prev[j] = 0.0
            f_prev[j] = 0.0

    # DP_table for visualization (if needed, like Perl's @DP)
    # dp_table_viz = []
    # if not is_reversed_for_hirschberg and m > 0:
    #     dp_table_viz.append(list(s_prev))

    # Fill rows
    for i in range(1, m + 1):
        # Initialize column 0 for current row i
        if MODE == 'global':
            s_cur[0] = -(GAP_OPEN + (i - 1) * GAP_EXTEND)
            e_cur[0] = NEG_INF  # E[i,0] cannot exist (gap in A before A starts)
            f_cur[0] = s_cur[0]  # F[i,0] ends with gap in B
        elif MODE == 'lcs':
            s_cur[0] = 0.0
            e_cur[0] = NEG_INF
            f_cur[0] = NEG_INF
        else:  # local
            s_cur[0] = 0.0
            e_cur[0] = 0.0
            f_cur[0] = 0.0

        for j in range(1, n + 1):
            char_a = seq_a_list[i - 1]
            char_b = seq_b_list[j - 1]

            sub_score = score_matrix_chars(char_a, char_b)
            if MODE == 'lcs':
                sub_score = 1 if char_a == char_b else 0

            # M_val: Score from diagonal precursor states S[i-1,j-1], E[i-1,j-1], F[i-1,j-1]
            diag_pred_max = max(s_prev[j - 1], e_prev[j - 1], f_prev[j - 1])
            m_val = diag_pred_max + sub_score

            # E_val: Score from left S[i,j-1] or E[i,j-1] (gap in seq_a_list, char_b vs '-')
            e_val_open = s_cur[j - 1] - GAP_OPEN
            e_val_ext = e_cur[j - 1] - GAP_EXTEND
            e_val = max(e_val_open, e_val_ext)

            # F_val: Score from top S[i-1,j] or F[i-1,j] (gap in seq_b_list, char_a vs '-')
            f_val_open = s_prev[j] - GAP_OPEN
            f_val_ext = f_prev[j] - GAP_EXTEND
            f_val = max(f_val_open, f_val_ext)

            if MODE == 'local':
                m_val = max(0, m_val)
                e_val = max(0, e_val)
                f_val = max(0, f_val)

            s_val = max(m_val, e_val, f_val)
            if MODE == 'local':  # S_val also max(0,...)
                s_val = max(0, s_val)

            s_cur[j] = s_val
            e_cur[j] = e_val
            f_cur[j] = f_val

        # if not is_reversed_for_hirschberg and m > 0:
        #     dp_table_viz.append(list(s_cur))

        s_prev, e_prev, f_prev = list(s_cur), list(e_cur), list(f_cur)  # Copy current to previous for next iteration

        # Progress bar (only if not in a deeply nested serial call within a fork)
        # The IN_CHILD check might need refinement for Python's multiprocessing
        if not IN_CHILD or IN_CHILD == 1:  # Simplistic check
            if not is_reversed_for_hirschberg:  # only for forward pass of hirschberg, or non-hirschberg
                progress_bar(i, m)

    return s_cur  # Return the last score row (S scores)


# --- Multiprocessing Wrapper (Python's maybe_fork) ---
# Global pool, initialized in main
_pool = None


def init_worker(gap_open_init, gap_extend_init, mode_init, matrix_init, neg_inf_init, in_child_init_val,
                fork_results_proxy):
    """Initializer for worker processes in the pool."""
    global GAP_OPEN, GAP_EXTEND, MODE, MATRIX, NEG_INF, IN_CHILD, FORK_RESULTS
    GAP_OPEN = gap_open_init
    GAP_EXTEND = gap_extend_init
    MODE = mode_init
    MATRIX = matrix_init
    NEG_INF = neg_inf_init
    IN_CHILD = in_child_init_val  # Workers are "children"
    FORK_RESULTS = fork_results_proxy  # This won't work as intended for direct result passing like Perl hash.
    # Results will be returned by apply_async.


def maybe_fork_python(task_id, function_ref, *args):
    """
    Executes function_ref(*args) either in a separate process or serially.
    Returns the result of the function.
    task_id is for potential debugging or specific result tracking if needed.
    """
    global IN_CHILD, _pool, FORK_RESULTS

    # Check if the task is too small or if we are already in a child/worker
    # or if MAX_PROCS is low.
    # For problem size, we'd need to inspect args, e.g. sequence lengths.
    # Assuming args[0] and args[1] are sequences or lists for simplicity here.
    # This heuristic needs to be adapted based on what function_ref is.
    current_problem_size = 1
    if function_ref == gotoh_linear_space_scores_wrapper or function_ref == small_align_wrapper:
        # args for gotoh_linear_space_scores_wrapper: (seq_a_list, seq_b_list, is_reversed)
        # args for small_align_wrapper: (seq_a_list, seq_b_list)
        seq_a_arg = args[0][0]  # The first argument to the *wrapper* is a tuple of args for the target
        seq_b_arg = args[0][1]
        current_problem_size = len(seq_a_arg) * len(seq_b_arg)

    if IN_CHILD or MAX_PROCS <= 1 or (_pool is None) or \
            (function_ref in [gotoh_linear_space_scores_wrapper,
                              small_align_wrapper] and current_problem_size < MIN_SIZE_FOR_FORK):
        # Run serially
        original_in_child = IN_CHILD
        IN_CHILD = True  # Mark that we are in a "child-like" execution path
        result = function_ref(*args)
        IN_CHILD = original_in_child
        return result

    # Run in parallel using the global pool
    # The FORK_RESULTS mechanism from Perl (shared hash updated by callback)
    # is tricky with multiprocessing.Pool. It's more common to get results directly.
    # We use apply_async and get the result.
    # For Hirschberg, we need two results (L and R) which could be one task returning a tuple,
    # or two separate async tasks.

    # If function_ref is designed to be called like hirschberg_L_R_task:
    if task_id.startswith("HB_"):  # Special handling for Hirschberg's combined L/R task
        async_result = _pool.apply_async(hirschberg_L_R_task_wrapper, args=args)
    else:  # For NWScore or small_align tasks
        async_result = _pool.apply_async(function_ref, args=args)

    try:
        # Add a timeout to prevent indefinite blocking
        # This timeout needs to be generous or configurable
        timeout_seconds = 3600  # 1 hour, example
        result = async_result.get(timeout=timeout_seconds)
        return result
    except Exception as e:
        print(f"Error in forked task {task_id}: {e}", file=sys.stderr)
        # Fallback or re-raise. For Hirschberg, a fallback to serial is in place.
        return None  # Indicate failure


# --- Wrappers for multiprocessing ---
# Functions passed to Pool.apply_async must be picklable (top-level).
# These wrappers take a tuple of arguments.

def gotoh_linear_space_scores_wrapper(args_tuple):
    # args_tuple contains (seq_a_list, seq_b_list, is_reversed)
    return gotoh_linear_space_scores(*args_tuple)


def small_align_wrapper(args_tuple):
    # args_tuple contains (seq_a_list, seq_b_list)
    return small_align_core(*args_tuple)  # small_align_core is the pure computation part


def hirschberg_L_R_task_wrapper(args_tuple):
    # args_tuple contains (seq_a_list_L, seq_b_list_L, is_reversed_L,
    #                      seq_a_list_R, seq_b_list_R, is_reversed_R)
    # This is specific to how hirschberg calls it.
    seq_a_L, seq_b_L, rev_L, seq_a_R, seq_b_R, rev_R = args_tuple

    # These are direct calls because they are inside a worker already.
    # We need to ensure IN_CHILD is True for them if they are not to fork further.
    # The worker init sets IN_CHILD=True.
    l_scores = gotoh_linear_space_scores(seq_a_L, seq_b_L, rev_L)
    r_scores = gotoh_linear_space_scores(seq_a_R, seq_b_R, rev_R)
    return l_scores, r_scores


# --- Alignment Algorithms ---
def nw_score_facade(seq_a_list, seq_b_list, is_reversed):
    """Facade for gotoh_linear_space_scores for use in Hirschberg, uses maybe_fork."""
    task_id = f"NW_{is_reversed}_{len(seq_a_list)}x{len(seq_b_list)}"
    # Pass arguments as a tuple to the wrapper
    return maybe_fork_python(task_id, gotoh_linear_space_scores_wrapper, (seq_a_list, seq_b_list, is_reversed))


def small_align_facade(seq_a_list, seq_b_list):
    """Facade for small_align (full DP with traceback), uses maybe_fork."""
    task_id = f"SA_{len(seq_a_list)}x{len(seq_b_list)}"
    return maybe_fork_python(task_id, small_align_wrapper, (seq_a_list, seq_b_list))


def small_align_core(seq_a_list, seq_b_list):
    """
    Corresponds to Perl's small_align core logic.
    Full DP with traceback for small sequences or base cases.
    """
    global IN_CHILD
    m = len(seq_a_list)
    n = len(seq_b_list)

    s_mat = [[0.0] * (n + 1) for _ in range(m + 1)]
    e_mat = [[NEG_INF] * (n + 1) for _ in range(m + 1)]  # Gap in A
    f_mat = [[NEG_INF] * (n + 1) for _ in range(m + 1)]  # Gap in B
    p_mat = [[''] * (n + 1) for _ in range(m + 1)]  # Pointer matrix

    s_mat[0][0] = 0.0
    p_mat[0][0] = 'X'  # Stopper

    # Initialize first row (i=0)
    for j in range(1, n + 1):
        if MODE == 'global' or MODE == 'lcs':
            val = -(GAP_OPEN + (j - 1) * GAP_EXTEND)
            s_mat[0][j] = val
            e_mat[0][j] = val
            f_mat[0][j] = NEG_INF
            p_mat[0][j] = 'E'  # Gap in A (horizontal)
        else:  # local
            s_mat[0][j], e_mat[0][j], f_mat[0][j] = 0.0, 0.0, 0.0
            p_mat[0][j] = 'X'

    # Initialize first column (j=0)
    for i in range(1, m + 1):
        if MODE == 'global' or MODE == 'lcs':
            val = -(GAP_OPEN + (i - 1) * GAP_EXTEND)
            s_mat[i][0] = val
            e_mat[i][0] = NEG_INF
            f_mat[i][0] = val
            p_mat[i][0] = 'F'  # Gap in B (vertical)
        else:  # local
            s_mat[i][0], e_mat[i][0], f_mat[i][0] = 0.0, 0.0, 0.0
            p_mat[i][0] = 'X'

    # Fill matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            char_a = seq_a_list[i - 1]
            char_b = seq_b_list[j - 1]
            sub = score_matrix_chars(char_a, char_b)
            if MODE == 'lcs':
                sub = 1 if char_a == char_b else 0

            # M_val comes from max of S,E,F of [i-1,j-1]
            diag_pred = max(s_mat[i - 1][j - 1], e_mat[i - 1][j - 1], f_mat[i - 1][j - 1])
            m_val_cand = diag_pred + sub

            # E_val (gap in A) from S[i,j-1] or E[i,j-1]
            e_val_cand = max(s_mat[i][j - 1] - GAP_OPEN, e_mat[i][j - 1] - GAP_EXTEND)

            # F_val (gap in B) from S[i-1,j] or F[i-1,j]
            f_val_cand = max(s_mat[i - 1][j] - GAP_OPEN, f_mat[i - 1][j] - GAP_EXTEND)

            if MODE == 'local':
                m_val_cand = max(0, m_val_cand)
                e_val_cand = max(0, e_val_cand)
                f_val_cand = max(0, f_val_cand)

            e_mat[i][j] = e_val_cand
            f_mat[i][j] = f_val_cand

            best_score = m_val_cand
            ptr = 'M'

            if e_val_cand > best_score:  # Python careful with >= vs > for pointer preference
                best_score = e_val_cand
                ptr = 'E'

            if f_val_cand > best_score:
                best_score = f_val_cand
                ptr = 'F'

            if MODE == 'local' and 0 > best_score:  # If all paths negative, local starts new
                s_mat[i][j] = 0
                p_mat[i][j] = 'X'
                e_mat[i][j] = max(0, e_mat[i][j])  # Ensure E,F also non-negative if S is 0
                f_mat[i][j] = max(0, f_mat[i][j])
            else:
                s_mat[i][j] = best_score
                p_mat[i][j] = ptr

        if not IN_CHILD or IN_CHILD == 1:
            progress_bar(i, m)

    # Traceback
    align_a, align_b = "", ""
    cur_i, cur_j = m, n

    if MODE == 'local':
        max_s = NEG_INF  # For local, find max score in S matrix
        # For local, an empty alignment has score 0. So init max_s to 0 or NEG_INF
        # If all scores are 0, (0,0) is fine. If positive score exists, find it.
        max_s = 0.0
        cur_i, cur_j = 0, 0  # Default if no positive alignment
        for r_idx in range(m + 1):
            for c_idx in range(n + 1):
                if s_mat[r_idx][c_idx] > max_s:
                    max_s = s_mat[r_idx][c_idx]
                    cur_i, cur_j = r_idx, c_idx

    while cur_i > 0 or cur_j > 0:
        if p_mat[cur_i][cur_j] == 'X':
            break
        if p_mat[cur_i][cur_j] == 'M':
            align_a = seq_a_list[cur_i - 1] + align_a
            align_b = seq_b_list[cur_j - 1] + align_b
            cur_i -= 1
            cur_j -= 1
        elif p_mat[cur_i][cur_j] == 'E':  # Gap in A
            align_a = '-' + align_a
            align_b = seq_b_list[cur_j - 1] + align_b
            cur_j -= 1
        elif p_mat[cur_i][cur_j] == 'F':  # Gap in B
            align_a = seq_a_list[cur_i - 1] + align_a
            align_b = '-' + align_b
            cur_i -= 1
        else:  # Should not happen with 'X' as stopper
            break

    return align_a, align_b


def hirschberg(seq_a_list, seq_b_list):
    """Hirschberg's divide-and-conquer alignment algorithm."""
    global IN_CHILD
    m = len(seq_a_list)
    n = len(seq_b_list)

    if m == 0:
        return '-' * n, "".join(seq_b_list)
    if n == 0:
        return "".join(seq_a_list), '-' * m
    if m == 1 or n == 1:
        # For base case, call the core logic directly, not the facade,
        # to avoid potential re-forking if not handled by IN_CHILD or MIN_SIZE_FOR_FORK in facade.
        # Or, ensure facade's serial path is robust. For now, direct call for base.
        original_in_child = IN_CHILD
        IN_CHILD = True  # Ensure serial execution for this small part
        res = small_align_core(seq_a_list, seq_b_list)
        IN_CHILD = original_in_child
        return res

    i_mid = m // 2
    task_id = f"HB_{m}x{n}"

    # Arguments for the combined L/R task
    # L: A_top vs B. A_top = seq_a_list[0:i_mid]
    # R: A_bottom_rev vs B_rev. A_bottom = seq_a_list[i_mid:m]
    #    A_bottom_rev = list(reversed(seq_a_list[i_mid:m]))
    #    B_rev = list(reversed(seq_b_list))
    # In Perl, B_rev was implicit in _NWScore_serial(..., rev=1).
    # Here, gotoh_linear_space_scores takes B as is, and if rev=1, reverses both A and B.
    # So for R_scores, we need to pass A_bottom_rev and B_rev to gotoh_linear_space_scores
    # and set its is_reversed_for_hirschberg to True (or handle reversal inside).

    # For R_scores(A_bottom, B, rev=1):
    #   A_bottom_rev = list(reversed(seq_a_list[i_mid:m]))
    #   B_rev for R_scores = list(reversed(seq_b_list))
    #   R_scores_raw = gotoh_linear_space_scores(A_bottom_rev, B_rev, True)
    #   Then R_scores_final[j] is R_scores_raw[n-j] (due to B being reversed).

    # The Perl script's _NWScore_serial with rev=1 reverses A and B.
    # So, L_scores uses A[0..i_mid-1] and B, not reversed.
    # R_scores uses A[i_mid..m-1] (reversed) and B (reversed).

    seq_a_L_part = seq_a_list[0:i_mid]
    seq_a_R_part_for_rev = seq_a_list[i_mid:m]  # This will be reversed by gotoh if rev=1

    # Call a task that computes L and R scores, possibly in parallel
    # The task hirschberg_L_R_task_wrapper does:
    #   l_scores = gotoh_linear_space_scores(seq_a_L_part, seq_b_list, False)
    #   r_scores = gotoh_linear_space_scores(seq_a_R_part_for_rev, seq_b_list, True)
    #                                       (where gotoh with True reverses both inputs)

    fork_res = maybe_fork_python(
        task_id,
        hirschberg_L_R_task_wrapper,
        (seq_a_L_part, seq_b_list, False,  # Args for L
         seq_a_R_part_for_rev, seq_b_list, True)  # Args for R
    )

    l_score_row, r_score_row = None, None
    if fork_res and isinstance(fork_res, tuple) and len(fork_res) == 2 and \
            isinstance(fork_res[0], list) and isinstance(fork_res[1], list):
        l_score_row, r_score_row = fork_res
    else:
        # Fallback to serial if fork failed or returned bad data
        if MAX_PROCS > 1 and _pool is not None:  # Only warn if parallel was attempted
            print(f"Fork for Hirschberg {task_id} failed/returned invalid, running L/R serially.", file=sys.stderr)

        original_in_child = IN_CHILD
        IN_CHILD = True  # Mark as "in child" to ensure serial execution of gotoh
        l_score_row = gotoh_linear_space_scores(seq_a_L_part, seq_b_list, False)
        # For r_score_row, gotoh_linear_space_scores with rev=True will reverse its inputs.
        r_score_row = gotoh_linear_space_scores(seq_a_R_part_for_rev, seq_b_list, True)
        IN_CHILD = original_in_child

    if not l_score_row or not r_score_row:
        raise Exception(f"Failed to compute score arrays for Hirschberg split {task_id}")

    best_score_sum = NEG_INF
    j_mid = 0
    for j in range(n + 1):
        # l_score_row is S(A_top, B_left_j)
        # r_score_row is S(A_bottom_rev, B_rev_right_of_j_split)
        #   which is S(A_bottom, B_right_of_j_split)
        #   B_right_of_j_split has length (n-j). Its reverse is used.
        #   So r_score_row[n-j] is what we need.
        score_l = l_score_row[j]
        score_r = r_score_row[n - j]  # B was reversed for r_score_row calculation.
        current_sum = score_l + score_r
        if current_sum > best_score_sum:
            best_score_sum = current_sum
            j_mid = j

    # Recursive calls
    # These will use their own maybe_fork logic
    a1, b1 = hirschberg(seq_a_list[0:i_mid], seq_b_list[0:j_mid])
    a2, b2 = hirschberg(seq_a_list[i_mid:m], seq_b_list[j_mid:n])

    return a1 + a2, b1 + b2


# --- Output and Utility Functions ---
def extract_lcs_from_alignment(aln_a, aln_b):
    """Extracts longest common subsequence runs from aligned sequences."""
    runs = []
    current_run = ""
    for char_a, char_b in zip(aln_a, aln_b):
        if char_a == char_b and char_a != '-':
            current_run += char_a
        else:
            if current_run:
                runs.append(current_run)
            current_run = ""
    if current_run:
        runs.append(current_run)
    return runs


def build_full_dp_matrix(seq_a_list, seq_b_list):
    """Builds the full S-score DP matrix (Python version of build_full_DP)."""
    m = len(seq_a_list)
    n = len(seq_b_list)

    full_dp_S = [[0.0] * (n + 1) for _ in range(m + 1)]

    s_prev = [0.0] * (n + 1)
    e_prev = [NEG_INF] * (n + 1)
    f_prev = [NEG_INF] * (n + 1)

    s_cur = [0.0] * (n + 1)
    e_cur = [0.0] * (n + 1)
    f_cur = [0.0] * (n + 1)

    # Init row 0 for s_prev, e_prev, f_prev (same as in gotoh_linear_space_scores)
    s_prev[0] = 0.0
    e_prev[0], f_prev[0] = NEG_INF, NEG_INF
    for j in range(1, n + 1):
        if MODE == 'global':
            s_prev[j] = -(GAP_OPEN + (j - 1) * GAP_EXTEND)
            e_prev[j], f_prev[j] = s_prev[j], NEG_INF
        elif MODE == 'lcs':
            s_prev[j], e_prev[j], f_prev[j] = 0.0, NEG_INF, NEG_INF
        else:  # local
            s_prev[j], e_prev[j], f_prev[j] = 0.0, 0.0, 0.0
    full_dp_S[0] = list(s_prev)

    # Fill rows
    for i in range(1, m + 1):
        # Init col 0 for s_cur, e_cur, f_cur
        if MODE == 'global':
            s_cur[0] = -(GAP_OPEN + (i - 1) * GAP_EXTEND)
            e_cur[0], f_cur[0] = NEG_INF, s_cur[0]
        elif MODE == 'lcs':
            s_cur[0], e_cur[0], f_cur[0] = 0.0, NEG_INF, NEG_INF
        else:  # local
            s_cur[0], e_cur[0], f_cur[0] = 0.0, 0.0, 0.0

        for j in range(1, n + 1):
            char_a, char_b = seq_a_list[i - 1], seq_b_list[j - 1]
            sub = score_matrix_chars(char_a, char_b)
            if MODE == 'lcs': sub = 1 if char_a == char_b else 0

            diag_pred = max(s_prev[j - 1], e_prev[j - 1], f_prev[j - 1])
            m_val = diag_pred + sub
            e_val = max(s_cur[j - 1] - GAP_OPEN, e_cur[j - 1] - GAP_EXTEND)
            f_val = max(s_prev[j] - GAP_OPEN, f_prev[j] - GAP_EXTEND)

            if MODE == 'local':
                m_val, e_val, f_val = max(0, m_val), max(0, e_val), max(0, f_val)

            s_val = max(m_val, e_val, f_val)
            if MODE == 'local': s_val = max(0, s_val)

            s_cur[j], e_cur[j], f_cur[j] = s_val, e_val, f_val

        full_dp_S[i] = list(s_cur)
        s_prev, e_prev, f_prev = list(s_cur), list(e_cur), list(f_cur)
        # No progress bar in this one for simplicity, assumed to be called once.

    return full_dp_S


def print_colored_alignment(aln_a, aln_b):
    """Prints alignment with colors."""
    len_aln = len(aln_a)
    if len_aln == 0: return

    pos_a, pos_b = 0, 0  # 0-based ungapped positions

    for offset in range(0, len_aln, 80):
        block_a = aln_a[offset: offset + 80]
        block_b = aln_b[offset: offset + 80]
        blen = len(block_a)

        # Seq A line
        sys.stdout.write(f"{pos_a + 1:6d} ")
        current_line_pos_a = pos_a
        for k in range(blen):
            ca, cb = block_a[k], block_b[k]
            color = COL_RED if ca == '-' or cb == '-' else \
                COL_GREEN if ca == cb else COL_CYAN
            sys.stdout.write(f"{color}{ca}{COL_RESET}")
            if ca != '-': current_line_pos_a += 1
        sys.stdout.write(f" {current_line_pos_a}\n")

        # Seq B line
        sys.stdout.write(f"{pos_b + 1:6d} ")
        current_line_pos_b = pos_b
        for k in range(blen):
            ca, cb = block_a[k], block_b[k]
            color = COL_RED if ca == '-' or cb == '-' else \
                COL_GREEN if ca == cb else COL_CYAN
            sys.stdout.write(f"{color}{cb}{COL_RESET}")
            if cb != '-': current_line_pos_b += 1
        sys.stdout.write(f" {current_line_pos_b}\n\n")

        pos_a = current_line_pos_a  # Update overall ungapped positions
        pos_b = current_line_pos_b
    sys.stdout.flush()


def progress_bar(done, total, width=40):
    """Prints a progress bar to STDERR."""
    if total == 0: return
    pct = done / total
    pct = min(1, max(0, pct))
    filled = int(pct * width)
    bar = "[" + "#" * filled + " " * (width - filled) + "]"
    sys.stderr.write(f"\r{bar} {pct * 100:3.0f}%")
    if done >= total:
        sys.stderr.write("\n")
    sys.stderr.flush()


# --- Main Application Logic ---
def main():
    global GAP_OPEN, GAP_EXTEND, MODE, MATRIX, NEG_INF, MAX_PROCS, _pool, IN_CHILD
    # _fork_results_manager, FORK_RESULTS # Not used in current Pythonic parallel approach

    parser = argparse.ArgumentParser(description="Pairwise sequence alignment tool (Python port).")
    parser.add_argument("seq_a_file", help="Input FASTA file for sequence A")
    parser.add_argument("seq_b_file", help="Input FASTA file for sequence B")
    parser.add_argument("--mode", choices=['global', 'local', 'lcs'], default='global',
                        help="Alignment mode (default: global)")
    parser.add_argument("--matrix", required=True, help="Scoring matrix file")
    parser.add_argument("--gapopen", type=int, default=10, help="Gap open penalty (positive value)")
    parser.add_argument("--gapext", type=int, default=1, help="Gap extend penalty (positive value)")
    parser.add_argument("--out", default="align", help="Output prefix for files")
    parser.add_argument("--maxprocs", type=int, default=MAX_PROCS, help="Max processes for parallel parts")

    args = parser.parse_args()

    MODE = args.mode
    GAP_OPEN = args.gapopen  # Store as positive, apply as negative in algorithm
    GAP_EXTEND = args.gapext
    MAX_PROCS = args.maxprocs

    try:
        MATRIX = read_matrix_file(args.matrix)
    except Exception as e:
        print(f"Error reading matrix file '{args.matrix}': {e}", file=sys.stderr)
        sys.exit(1)

    seq_a_str = slurp_seq(args.seq_a_file)
    seq_b_str = slurp_seq(args.seq_b_file)
    seq_a_list = list(seq_a_str)
    seq_b_list = list(seq_b_str)

    if MAX_PROCS > 1:
        # Initialize multiprocessing pool
        # The initializer passes necessary global-like variables to worker processes.
        # However, for complex shared state like FORK_RESULTS hash, it's better if tasks return results.
        init_args = (GAP_OPEN, GAP_EXTEND, MODE, MATRIX, NEG_INF, True, None)  # Pass True for IN_CHILD in workers
        _pool = Pool(processes=MAX_PROCS, initializer=init_worker, initargs=init_args)
    else:
        _pool = None  # Ensure serial execution if MAX_PROCS <=1

    aligned_s1, aligned_s2 = "", ""
    if MODE == 'lcs':
        # LCS uses full DP (small_align_core), no Hirschberg typically
        # Call facade to respect potential forking for very large LCS base cases if MIN_SIZE_FOR_FORK is met
        # aligned_s1, aligned_s2 = small_align_facade(seq_a_list, seq_b_list)
        # Or directly if LCS is always non-Hirschberg in Perl script:
        IN_CHILD = True  # Simulate being in a non-forking path for this direct call
        aligned_s1, aligned_s2 = small_align_core(seq_a_list, seq_b_list)
        IN_CHILD = False
    else:  # global or local
        aligned_s1, aligned_s2 = hirschberg(seq_a_list, seq_b_list)

    if _pool:
        _pool.close()
        _pool.join()

    print("\n")  # After progress bars
    print_colored_alignment(aligned_s1, aligned_s2)

    # --- Output files ---
    with open(f"{args.out}.A.fa", "w") as f:
        f.write(f">A_aligned_{MODE}\n{aligned_s1}\n")
    with open(f"{args.out}.B.fa", "w") as f:
        f.write(f">B_aligned_{MODE}\n{aligned_s2}\n")

    print(f"Aligned sequences saved to {args.out}.A.fa and {args.out}.B.fa", file=sys.stderr)

    dp_full = build_full_dp_matrix(seq_a_list, seq_b_list)

    # Write TSV
    try:
        with open(f"{args.out}.matrix.tsv", "w") as f_tsv:
            f_tsv.write("\t" + "\t".join(seq_b_list) + "\n")
            for i, row_scores in enumerate(dp_full):
                label = seq_a_list[i - 1] if i > 0 else ''
                f_tsv.write(label + "\t" + "\t".join(map(str, row_scores)) + "\n")
        print(f"Full DP matrix saved to {args.out}.matrix.tsv", file=sys.stderr)
    except Exception as e:
        print(f"Error writing TSV matrix: {e}", file=sys.stderr)

    # Binary DP matrix (omitted for brevity, similar packing logic needed if required)

    if MODE == 'lcs':
        lcs_runs = extract_lcs_from_alignment(aligned_s1, aligned_s2)
        if lcs_runs:
            longest_lcs = max(lcs_runs, key=len)
            with open(f"{args.out}.lcs.fa", "w") as f_lcs:
                f_lcs.write(f">{args.out}_LCS\n{longest_lcs}\n")
            print(f"Longest LCS saved to {args.out}.lcs.fa", file=sys.stderr)
        else:
            with open(f"{args.out}.lcs.fa", "w") as f_lcs:  # Create empty if none
                f_lcs.write(f">{args.out}_LCS\n\n")
            print(f"No common subsequence found for {args.out}.lcs.fa", file=sys.stderr)

    final_score = "N/A"
    if dp_full:
        m_dp, n_dp = len(dp_full) - 1, len(dp_full[0]) - 1
        if MODE == 'local':
            final_score = 0.0
            for r in dp_full:
                for cell_score in r:
                    if cell_score > final_score:
                        final_score = cell_score
        else:  # global or lcs
            final_score = dp_full[m_dp][n_dp]

    print(f"\nScore   :  {final_score}")
    print(f"Output  :  {args.out}.A.fa, {args.out}.B.fa")
    print(f"Completed {MODE} alignment.")


if __name__ == "__main__":
    main()