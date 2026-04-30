import math
from src.hmm import DELETE, MATCH

def viterbi(seq, emissions, transitions):
    states = [k for k in transitions if k not in ("Start", "End")]

    m = []
    path = {} # We will actively carry the winning paths here

    # --- COLUMN 0 ---
    first_col = {}
    for s in states:
        prob = transitions["Start"].get(s, 0)
        if prob > 0 and s[0] != DELETE and seq[0] in emissions[s]:
            first_col[s] = math.log(prob) + math.log(emissions[s][seq[0]])
        else:
            first_col[s] = -1e9
            
        # Initialize the path history for the first column
        path[s] = [s] 
        
    m.append(first_col)

    # --- COLUMNS 1 to N ---
    for i in range(1, len(seq)):
        char = seq[i]
        curr_col = {}
        newpath = {} # <--- Tracks the paths for the CURRENT column
        
        # HELPER: Now returns (score, winning_previous_state)
        def get_best_jump(target, possible_prevs, prev_matrix, is_delete=False):
            candidates = []
            for p in possible_prevs:
                if p in prev_matrix and (jump := transitions[p].get(target, 0)) > 0:
                    score = prev_matrix[p] + math.log(jump)
                    if not is_delete and char in emissions[target]:
                        score += math.log(emissions[target][char])
                    candidates.append((score, p))
            return max(candidates, key=lambda x: x[0]) if candidates else (-1e9, None)

        # PASS 1: Calculate M and I (looking at i-1)
        for s in states:
            stype, num = s[0], int(s[1:])
            if stype == MATCH:
                score, best_p = get_best_jump(s, [f"M{num-1}", f"I{num-1}", f"D{num-1}"], m[i-1])
            elif stype == "I":
                score, best_p = get_best_jump(s, [f"M{num}", s, f"D{num}"], m[i-1])
            else:
                continue
                
            curr_col[s] = score
            
            # THE FIX: Safely get the path, or default to an empty list if None/missing
            newpath[s] = path.get(best_p, []) + [s]

        # PASS 2: Calculate D (looking at the column we JUST built above)
        for s in states:
            if s[0] == DELETE:
                num = int(s[1:])
                score, best_p = get_best_jump(s, [f"M{num-1}", f"D{num-1}"], curr_col, is_delete=True)
                
                curr_col[s] = score
                
                # THE FIX: Safely get from newpath since D states look at the current column
                newpath[s] = newpath.get(best_p, []) + [s]

        m.append(curr_col)
        
        # Overwrite the old paths to save memory, just like GfG
        path = newpath 

    # Because we carried the paths with us, we don't need a traceback loop!
    best_final_score = -1e9
    best_last_state = None

    for state, score in m[-1].items():
        jump_to_end = transitions[state].get("End", 0)
        
        if jump_to_end > 0:
            final_score = score + math.log(jump_to_end)
            if final_score > best_final_score:
                best_final_score = final_score
                best_last_state = state

    # Just return the winning path that we carried all the way to the end
    return path.get(best_last_state, [])
