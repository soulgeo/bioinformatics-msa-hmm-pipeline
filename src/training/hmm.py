import math
from collections import defaultdict

from src.generation.generate import ALPHABET
from src.settings import DELETE, INSERT, MATCH


class HMM:
    def __init__(self, num_matches, emissions, transitions):
        self.num_matches = num_matches
        self.emissions = emissions
        self.transitions = transitions

    @classmethod
    def from_msa(cls, msa):
        states_seq = cls._create_states_sequence(msa)
        num_matches = states_seq.count(MATCH)
        emissions = cls._calculate_emissions_matrix(msa, states_seq)
        transitions = cls._calculate_transitions_matrix(msa, states_seq)
        return cls(num_matches, emissions, transitions)

    def decode(self, gene):
        """Finds the most likely path for a gene (Viterbi algorithm)."""
        states = [k for k in self.transitions if k not in ("Start", "End")]
        m = []
        path = {}

        # --- COLUMN 0 ---
        first_col = {}
        for state in states:
            s_type = state[0]
            prob = self.transitions["Start"].get(state, 0)
            if (
                prob > 0
                and s_type != DELETE
                and gene[0] in self.emissions[state]
            ):
                first_col[state] = math.log(prob) + math.log(
                    self.emissions[state][gene[0]]
                )
            else:
                first_col[state] = -1e9
            path[state] = [state]
        m.append(first_col)

        # --- COLUMNS 1 to N ---
        for i in range(1, len(gene)):
            char = gene[i]
            curr_col = {}
            newpath = {}

            def get_best_jump(target, possible_prevs, prev_matrix, is_delete=False):
                candidates = []
                for p in possible_prevs:
                    if (
                        p in prev_matrix
                        and (jump := self.transitions[p].get(target, 0)) > 0
                    ):
                        score = prev_matrix[p] + math.log(jump)
                        if not is_delete and char in self.emissions[target]:
                            score += math.log(self.emissions[target][char])
                        candidates.append((score, p))
                return (
                    max(candidates, key=lambda x: x[0])
                    if candidates
                    else (-1e9, None)
                )

            # PASS 1: Calculate M and I
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

            # PASS 2: Calculate D
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
            path = newpath

        best_final_score = -1e9
        best_last_state = None
        for state, score in m[-1].items():
            jump_to_end = self.transitions[state].get("End", 0)
            if jump_to_end > 0:
                final_score = score + math.log(jump_to_end)
                if final_score > best_final_score:
                    best_final_score = final_score
                    best_last_state = state

        return path.get(best_last_state, [])

    def forward(self, gene):
        """Calculates the forward score for a gene (Forward algorithm)."""
        states = [k for k in self.transitions if k not in ("Start", "End")]
        m = []

        # --- COLUMN 0 ---
        first_col = {}
        for state in states:
            prob = self.transitions["Start"].get(state, 0)
            if (
                prob > 0
                and state[0] != DELETE
                and gene[0] in self.emissions[state]
            ):
                first_col[state] = math.log(prob) + math.log(
                    self.emissions[state][gene[0]]
                )
            else:
                first_col[state] = -1e9
        m.append(first_col)

        # --- COLUMNS 1 to N ---
        for i in range(1, len(gene)):
            char = gene[i]
            curr_col = {}

            def get_total_jump(target, possible_prevs, prev_matrix, is_delete=False):
                candidates = []
                for p in possible_prevs:
                    if (
                        p in prev_matrix
                        and (jump := self.transitions[p].get(target, 0)) > 0
                    ):
                        score = prev_matrix[p] + math.log(jump)
                        if not is_delete and char in self.emissions[target]:
                            score += math.log(self.emissions[target][char])
                        candidates.append(score)
                return self._log_sum_exp(candidates)

            # PASS 1: Calculate M and I
            for state in states:
                s_type, num = state[0], int(state[1:])
                if s_type == MATCH:
                    curr_col[state] = get_total_jump(
                        state, [f"M{num-1}", f"I{num-1}", f"D{num-1}"], m[i - 1]
                    )
                elif s_type == INSERT:
                    curr_col[state] = get_total_jump(
                        state, [f"M{num}", state, f"D{num}"], m[i - 1]
                    )

            # PASS 2: Calculate D
            for state in states:
                if state[0] == DELETE:
                    num = int(state[1:])
                    curr_col[state] = get_total_jump(
                        state, [f"M{num-1}", f"D{num-1}"], curr_col, is_delete=True
                    )
            m.append(curr_col)

        # --- TERMINATION ---
        final_candidates = []
        for state, score in m[-1].items():
            jump_to_end = self.transitions[state].get("End", 0)
            if jump_to_end > 0:
                final_candidates.append(score + math.log(jump_to_end))

        return self._log_sum_exp(final_candidates)

    def train(self, dataset):
        """Performs Viterbi training on a dataset."""
        paths = [self.decode(gene) for gene in dataset]
        self._retrain_emissions(dataset, paths)
        self._retrain_transitions(paths)

    def _retrain_emissions(self, dataset, paths):
        emissions_tally = defaultdict(lambda: defaultdict(float))
        for seq, path in zip(dataset, paths):
            char_index = 0
            for state in path:
                if state.startswith("D"):
                    continue
                char = seq[char_index]
                emissions_tally[state][char] += 1
                char_index += 1

        pseudocount = 0.01
        all_possible_states = [f"I{k}" for k in range(self.num_matches + 1)] + [
            f"M{k}" for k in range(1, self.num_matches + 1)
        ]

        final_emissions = {}
        for state in all_possible_states:
            counts = emissions_tally[state]
            total_chars = sum(counts.values())
            divisor = total_chars + (len(ALPHABET) * pseudocount)
            final_emissions[state] = {
                letter: (counts.get(letter, 0) + pseudocount) / divisor
                for letter in ALPHABET
            }
        self.emissions = final_emissions

    def _retrain_transitions(self, paths):
        transitions_tally = {}
        # Pre-initialize topology
        transitions_tally["Start"] = {"M1": 0.0, "D1": 0.0, "I0": 0.0}
        for k in range(self.num_matches + 1):
            transitions_tally[f"I{k}"] = {f"I{k}": 0.0}
            if k < self.num_matches:
                transitions_tally[f"I{k}"][f"M{k+1}"] = 0.0
                transitions_tally[f"I{k}"][f"D{k+1}"] = 0.0
            else:
                transitions_tally[f"I{k}"]["End"] = 0.0
            if k > 0:
                transitions_tally[f"M{k}"] = {f"I{k}": 0.0}
                transitions_tally[f"D{k}"] = {f"I{k}": 0.0}
                if k < self.num_matches:
                    transitions_tally[f"M{k}"][f"M{k+1}"] = 0.0
                    transitions_tally[f"M{k}"][f"D{k+1}"] = 0.0
                    transitions_tally[f"D{k}"][f"M{k+1}"] = 0.0
                    transitions_tally[f"D{k}"][f"D{k+1}"] = 0.0
                else:
                    transitions_tally[f"M{k}"]["End"] = 0.0
                    transitions_tally[f"D{k}"]["End"] = 0.0

        for path in paths:
            prev_state = "Start"
            for current_state in path:
                transitions_tally[prev_state][current_state] += 1
                prev_state = current_state
            transitions_tally[prev_state]["End"] += 1

        pseudocount = 0.01
        for state, counts in transitions_tally.items():
            total_transitions = sum(counts.values())
            valid_paths = len(counts)
            divisor = total_transitions + (valid_paths * pseudocount)
            for next_state in counts.keys():
                prob = (counts[next_state] + pseudocount) / divisor
                transitions_tally[state][next_state] = prob
        self.transitions = transitions_tally

    @staticmethod
    def _log_sum_exp(log_probs):
        if not log_probs:
            return -1e9
        max_p = max(log_probs)
        if max_p <= -1e9:
            return -1e9
        sum_exp = sum(math.exp(p - max_p) for p in log_probs if p > -1e9)
        return max_p + math.log(sum_exp)

    @staticmethod
    def _create_states_sequence(msa):
        states_seq = []
        cols = len(msa[0])
        for i in range(cols):
            gaps = sum(1 for seq in msa if seq[i] == "-")
            states_seq.append(MATCH if gaps < len(msa) // 2 else INSERT)
        return states_seq

    @staticmethod
    def _calculate_emissions_matrix(msa, states_seq):
        emissions = defaultdict(lambda: defaultdict(float))
        match_counter = 0
        for i in range(len(states_seq)):
            current_type = states_seq[i]
            if current_type == MATCH:
                match_counter += 1
                state_key = f"M{match_counter}"
            else:
                state_key = f"I{match_counter}"
            for seq in msa:
                letter = seq[i]
                if letter in ALPHABET:
                    emissions[state_key][letter] += 1

        pseudocount = 0.01
        num_matches = states_seq.count(MATCH)
        all_possible_states = [f"I{k}" for k in range(num_matches + 1)] + [
            f"M{k}" for k in range(1, num_matches + 1)
        ]
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

    @staticmethod
    def _calculate_transitions_matrix(msa, states_seq):
        transitions = {}
        num_matches = states_seq.count(MATCH)
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

        pseudocount = 0.01
        for state, counts in transitions.items():
            total_transitions = sum(counts.values())
            valid_paths = len(counts)
            divisor = total_transitions + (valid_paths * pseudocount)
            for next_state in counts.keys():
                prob = (counts[next_state] + pseudocount) / divisor
                transitions[state][next_state] = prob
        return transitions
