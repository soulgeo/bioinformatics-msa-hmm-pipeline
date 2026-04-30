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
    emissions = {}

    # Step 1: Pre-initialize ALL required states
    num_matches = states_seq.count(MATCH)

    for k in range(num_matches + 1):
        emissions[f"I{k}"] = {letter: 0 for letter in ALPHABET}
        if k > 0:
            emissions[f"M{k}"] = {letter: 0 for letter in ALPHABET}

    # Step 2: Calculate the tallies
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

    # Step 3: Convert tallies to probabilities
    # probabilities of 0 will break Viterbi algorithm, so add a small pseudocount
    pseudocount = 0.01

    for state, counts in emissions.items():
        total_chars = sum(counts.values())
        divisor = total_chars + (len(ALPHABET) * pseudocount)

        for letter in ALPHABET:
            prob = (counts[letter] + pseudocount) / divisor
            emissions[state][letter] = prob

    return emissions


def calculate_transitions_matrix(msa, states_seq):
    transitions = defaultdict(lambda: defaultdict(float))

    # Step 1: Calculate the tallies
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

            # Found a gap, go through DELETE state if in MATCH, ignore if in INSERT
            if current_type == MATCH:
                current_state = f"D{state_counter}"
                transitions[temp_state][current_state] += 1
                temp_state = current_state
                
        # Transition to an "End" state when the sequence finishes
        transitions[temp_state]["End"] += 1

    # Step 2: Convert tallies to probabilities
    pseudocount = 0.01

    for state, counts in transitions.items():
        total_transitions = sum(counts.values())
        
        # Only apply the pseudocount multiplier to the valid outgoing transitions!
        valid_paths = len(counts) 
        divisor = total_transitions + (valid_paths * pseudocount)

        for next_state in counts.keys():
            prob = (counts[next_state] + pseudocount) / divisor
            transitions[state][next_state] = prob

    return transitions
