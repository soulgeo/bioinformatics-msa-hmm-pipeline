from src.generate import ALPHABET

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


def calculate_transition_matrix(msa, states_seq):
    pass
