import math

from src.training.hmm import DELETE, INSERT, MATCH

def log_sum_exp(log_probs):
    # If the list is empty, the probability is 0 (which is -inf in log space)
    if not log_probs:
        return -1e9
        
    # Find the biggest number to factor out
    max_p = max(log_probs)
    
    # If even the max is our "negative infinity", return it
    if max_p <= -1e9:
        return -1e9
        
    # The Log-Sum-Exp formula
    sum_exp = sum(math.exp(p - max_p) for p in log_probs if p > -1e9)
    return max_p + math.log(sum_exp)


def viterbi(gene, emissions, transitions):
    states = [k for k in transitions if k not in ("Start", "End")]

    # We are storing the "states x chars" viterbi matrix as a list of dictionaries.
    # Each dict represents a column of the matrix. The column contains the probs
    # for each state, for seq at that index. All probs will be calculated as logs
    # to avoid precision errors.

    m = []
    path = {}  # We will actively carry the winning paths here

    # --- COLUMN 0 ---
    first_col = {}
    for state in states:
        s_type = state[0]
        prob = transitions["Start"].get(state, 0)
        if prob > 0 and s_type != DELETE and gene[0] in emissions[state]:
            first_col[state] = math.log(prob) + math.log(
                emissions[state][gene[0]]
            )
        else:
            first_col[state] = -1e9 # Acts as -infinity for impossible paths

        # Initialize the path history for the first column
        path[state] = [state]

    m.append(first_col)

    # --- COLUMNS 1 to N ---
    for i in range(1, len(gene)):
        char = gene[i]
        curr_col = {}
        newpath = {}  # Tracks the paths for the CURRENT column

        # HELPER: Returns (score, winning_previous_state)
        def get_best_jump(target, possible_prevs, prev_matrix, is_delete=False):
            candidates = []
            for p in possible_prevs:
                if (
                    p in prev_matrix
                    and (jump := transitions[p].get(target, 0)) > 0
                ):
                    score = prev_matrix[p] + math.log(jump)
                    if not is_delete and char in emissions[target]:
                        score += math.log(emissions[target][char])
                    candidates.append((score, p))
            return (
                max(candidates, key=lambda x: x[0])
                if candidates
                else (-1e9, None)
            )

        # PASS 1: Calculate M and I (looking at i-1)
        for state in states:
            s_type, num = state[0], int(state[1:])
            if s_type == MATCH:
                possible_prevs = [f"M{num-1}", f"I{num-1}", f"D{num-1}"]
                score, best_p = get_best_jump(state, possible_prevs, m[i - 1])
            elif s_type == INSERT:
                possible_prevs = [f"M{num}", state, f"D{num}"]
                score, best_p = get_best_jump(state, possible_prevs, m[i - 1])
            else:
                continue

            curr_col[state] = score

            newpath[state] = path.get(best_p, []) + [state]

        # PASS 2: Calculate D (looking at the column we just built above)
        for state in states:
            if state[0] == DELETE:
                num = int(state[1:])
                possible_prevs = [f"M{num-1}", f"D{num-1}"]
                score, best_p = get_best_jump(
                    state, possible_prevs, curr_col, is_delete=True
                )

                curr_col[state] = score
                newpath[state] = newpath.get(best_p, []) + [state]

        m.append(curr_col)

        # Overwrite the old paths to save memory
        path = newpath

    best_final_score = -1e9
    best_last_state = None

    for state, score in m[-1].items():
        jump_to_end = transitions[state].get("End", 0)

        if jump_to_end > 0:
            final_score = score + math.log(jump_to_end)
            if final_score > best_final_score:
                best_final_score = final_score
                best_last_state = state

    # Return the winning path that we carried all the way to the end
    return path.get(best_last_state, [])


def forward_algorithm(gene, emissions, transitions):
    states = [k for k in transitions if k not in ("Start", "End")]
    m = []

    # --- COLUMN 0 ---
    first_col = {}
    for state in states:
        prob = transitions["Start"].get(state, 0)
        if prob > 0 and state[0] != DELETE and gene[0] in emissions[state]:
            first_col[state] = math.log(prob) + math.log(emissions[state][gene[0]])
        else:
            first_col[state] = -1e9
    m.append(first_col)

    # --- COLUMNS 1 to N ---
    for i in range(1, len(gene)):
        char = gene[i]
        curr_col = {}

        # HELPER: Returns the total summed score instead of the max score
        def get_total_jump(target, possible_prevs, prev_matrix, is_delete=False):
            candidates = []
            for p in possible_prevs:
                if p in prev_matrix and (jump := transitions[p].get(target, 0)) > 0:
                    score = prev_matrix[p] + math.log(jump)
                    if not is_delete and char in emissions[target]:
                        score += math.log(emissions[target][char])
                    candidates.append(score)
            
            return log_sum_exp(candidates)

        # PASS 1: Calculate M and I
        for state in states:
            s_type, num = state[0], int(state[1:])
            if s_type == MATCH:
                curr_col[state] = get_total_jump(state, [f"M{num-1}", f"I{num-1}", f"D{num-1}"], m[i - 1])
            elif s_type == INSERT:
                curr_col[state] = get_total_jump(state, [f"M{num}", state, f"D{num}"], m[i - 1])

        # PASS 2: Calculate D
        for state in states:
            if state[0] == DELETE:
                num = int(state[1:])
                curr_col[state] = get_total_jump(state, [f"M{num-1}", f"D{num-1}"], curr_col, is_delete=True)

        m.append(curr_col)

    # --- TERMINATION (Get the final sequence score) ---
    final_candidates = []
    
    for state, score in m[-1].items():
        jump_to_end = transitions[state].get("End", 0)
        if jump_to_end > 0:
            final_candidates.append(score + math.log(jump_to_end))

    # The total sum of all paths to the End state is our final sequence score!
    return log_sum_exp(final_candidates)
