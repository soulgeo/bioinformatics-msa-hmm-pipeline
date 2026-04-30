from collections import defaultdict
from src.generate import ALPHABET

def retrain_emissions_matrix(dataset_b, paths, num_matches):
    emissions = defaultdict(lambda: defaultdict(float))

    # Step 1: Calculate the tallies using the unique paths
    for seq, path in zip(dataset_b, paths):
        char_index = 0
        
        for state in path:
            # Delete states are silent, skip to the next state without moving char_index
            if state.startswith("D"):
                continue
                
            # M and I states emit the current character
            char = seq[char_index]
            emissions[state][char] += 1
            char_index += 1  # Only advance the character when an emission happens!

    # Step 2: Convert to probabilities (Same strict logic as before)
    pseudocount = 0.01
    all_possible_states = [f"I{k}" for k in range(num_matches + 1)] + \
                          [f"M{k}" for k in range(1, num_matches + 1)]

    final_emissions = {}
    for state in all_possible_states:
        counts = emissions[state]
        total_chars = sum(counts.values())
        divisor = total_chars + (len(ALPHABET) * pseudocount)
        
        final_emissions[state] = {
            letter: (counts.get(letter, 0) + pseudocount) / divisor 
            for letter in ALPHABET
        }
            
    return final_emissions


def retrain_transitions_matrix(paths, num_matches):
    transitions = {}

    # Step 1: Pre-initialize the exact same strict HMM topology as before
    transitions["Start"] = {"M1": 0.0, "D1": 0.0, "I0": 0.0}
    for k in range(num_matches + 1):
        transitions[f"I{k}"] = {f"I{k}": 0.0}
        if k < num_matches:
            transitions[f"I{k}"][f"M{k+1}"] = 0.0
            transitions[f"I{k}"][f"D{k+1}"] = 0.0
        else:
            transitions[f"I{k}"]["End"] = 0.0

        if k > 0:
            transitions[f"M{k}"] = {f"I{k}": 0.0}
            transitions[f"D{k}"] = {f"I{k}": 0.0}
            if k < num_matches:
                transitions[f"M{k}"][f"M{k+1}"] = 0.0
                transitions[f"M{k}"][f"D{k+1}"] = 0.0
                transitions[f"D{k}"][f"M{k+1}"] = 0.0
                transitions[f"D{k}"][f"D{k+1}"] = 0.0
            else:
                transitions[f"M{k}"]["End"] = 0.0
                transitions[f"D{k}"]["End"] = 0.0

    # Step 2: Calculate the tallies by simply reading the paths
    for path in paths:
        prev_state = "Start"
        for current_state in path:
            transitions[prev_state][current_state] += 1
            prev_state = current_state
            
        # Don't forget the final jump to End!
        transitions[prev_state]["End"] += 1

    # Step 3: Convert to probabilities
    pseudocount = 0.01
    for state, counts in transitions.items():
        total_transitions = sum(counts.values())
        valid_paths = len(counts) 
        divisor = total_transitions + (valid_paths * pseudocount)

        for next_state in counts.keys():
            prob = (counts[next_state] + pseudocount) / divisor
            transitions[state][next_state] = prob

    return transitions
