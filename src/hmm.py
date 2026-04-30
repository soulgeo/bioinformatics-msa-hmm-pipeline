from src.generate import ALPHABET
from collections import defaultdict

MATCH = "M"
INSERT = "I"
DELETE = "D"


def create_states_sequence(msa):
    states_seq = []
    cols = len(msa[0])
    for i in range(cols):
        gaps = 0
        for seq in msa:
            gaps += 1 if seq[i] == "-" else 0

        states_seq.append(MATCH if gaps < len(msa) // 2 else INSERT)

    return states_seq


def calculate_emissions_matrix(msa, states_seq):
    emissions = defaultdict(lambda: defaultdict(float))

    # Step 1: Calculate the tallies
    match_counter = 0

    # Iterate over all columns and add letters directly to their pre-made buckets
    for i in range(len(states_seq)):
        current_type = states_seq[i]

        if current_type == MATCH:
            match_counter += 1
            state_key = f"M{match_counter}"
        else:
            state_key = f"I{match_counter}"

        # Tally all the letters in the column
        for seq in msa:
            letter = seq[i]
            if letter in ALPHABET:
                emissions[state_key][letter] += 1

    # Step 2: Convert tallies to probabilities
    # Probabilities of 0 will break Viterbi algorithm, so add a small pseudocount
    pseudocount = 0.01

    # Make sure all states exist in the final matrix, even unobserved ones.
    num_matches = states_seq.count(MATCH)
    all_possible_states = [f"I{k}" for k in range(num_matches + 1)] + \
                         [f"M{k}" for k in range(1, num_matches + 1)]

    final_emissions = {}
    for state in all_possible_states:
        counts = emissions[state]
        total_chars = sum(counts.values())
        divisor = total_chars + (len(ALPHABET) * pseudocount)
        
        final_emissions[state] = {
            letter: (counts[letter] + pseudocount) / divisor 
            for letter in ALPHABET
        }
            
    return final_emissions


def calculate_transitions_matrix(msa, states_seq):
    transitions = {}
    num_matches = states_seq.count(MATCH)

    # Step 1: Pre-initialize the strict HMM topology
    # 1a. Start state
    transitions["Start"] = {"M1": 0.0, "D1": 0.0, "I0": 0.0}

    # 1b. Build the connections for all I, M, and D states
    for k in range(num_matches + 1):
        
        # Outgoing from Insert states
        transitions[f"I{k}"] = {f"I{k}": 0.0} # Self-loop
        if k < num_matches:
            transitions[f"I{k}"][f"M{k+1}"] = 0.0
            transitions[f"I{k}"][f"D{k+1}"] = 0.0
        else:
            transitions[f"I{k}"]["End"] = 0.0

        # Outgoing from Match and Delete states
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

    # Step 2: Calculate the tallies (Your exact logic here)
    for seq in msa:
        temp_state = "Start" 
        state_counter = 0
        
        for i, char in enumerate(seq):
            current_type = states_seq[i]
            if current_type == MATCH:
                state_counter += 1
                
            if char != "-":
                current_state = f"{current_type}{state_counter}"
                transitions[temp_state][current_state] += 1
                temp_state = current_state
                continue

            if current_type == MATCH:
                current_state = f"D{state_counter}"
                transitions[temp_state][current_state] += 1
                temp_state = current_state
                
        transitions[temp_state]["End"] += 1

    # Step 3: Convert tallies to probabilities
    pseudocount = 0.01

    for state, counts in transitions.items():
        total_transitions = sum(counts.values())
        
        valid_paths = len(counts) 
        divisor = total_transitions + (valid_paths * pseudocount)

        for next_state in counts.keys():
            prob = (counts[next_state] + pseudocount) / divisor
            transitions[state][next_state] = prob

    return transitions
