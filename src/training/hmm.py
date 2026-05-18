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

        # We temporarily create an instance to use helper methods
        instance = cls(num_matches, {}, {})
        instance.emissions = instance._calculate_emissions_matrix(
            msa, states_seq
        )
        instance.transitions = instance._calculate_transitions_matrix(
            msa, states_seq
        )
        return instance

    def decode(self, gene):
        """Finds the most likely path for a gene (Viterbi algorithm)."""
        states = [k for k in self.transitions if k not in ("Start", "End")]
        m = []
        path = {}

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
            path[state] = [state]
        m.append(first_col)

        # --- COLUMNS 1 to N ---
        for i in range(1, len(gene)):
            char = gene[i]
            curr_col, newpath = {}, {}

            def get_best_jump(
                target, possible_prevs, prev_matrix, is_delete=False
            ):
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
                    prevs = [f"M{num-1}", f"I{num-1}", f"D{num-1}"]
                elif s_type == INSERT:
                    prevs = [f"M{num}", state, f"D{num}"]
                else:
                    continue

                score, best_p = get_best_jump(state, prevs, m[i - 1])
                curr_col[state], newpath[state] = score, path.get(
                    best_p, []
                ) + [state]

            # PASS 2: Calculate D
            for state in states:
                if state[0] == DELETE:
                    num = int(state[1:])
                    score, best_p = get_best_jump(
                        state,
                        [f"M{num-1}", f"D{num-1}"],
                        curr_col,
                        is_delete=True,
                    )
                    curr_col[state], newpath[state] = score, newpath.get(
                        best_p, []
                    ) + [state]

            m.append(curr_col)
            path = newpath

        # Find best end state
        best_score, last_state = -1e9, None
        for state, score in m[-1].items():
            if (jump := self.transitions[state].get("End", 0)) > 0:
                final_score = score + math.log(jump)
                if final_score > best_score:
                    best_score, last_state = final_score, state

        return path.get(last_state, [])

    def forward(self, gene):
        """Calculates the forward score for a gene (Forward algorithm)."""
        states = [k for k in self.transitions if k not in ("Start", "End")]
        m = []

        # --- COLUMN 0 ---
        first_col = {s: -1e9 for s in states}
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
        m.append(first_col)

        # --- COLUMNS 1 to N ---
        for i in range(1, len(gene)):
            char = gene[i]
            curr_col = {}

            def get_total_jump(
                target, possible_prevs, prev_matrix, is_delete=False
            ):
                scores = []
                for p in possible_prevs:
                    if (
                        p in prev_matrix
                        and (jump := self.transitions[p].get(target, 0)) > 0
                    ):
                        score = prev_matrix[p] + math.log(jump)
                        if not is_delete and char in self.emissions[target]:
                            score += math.log(self.emissions[target][char])
                        scores.append(score)
                return self._log_sum_exp(scores)

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
                        state,
                        [f"M{num-1}", f"D{num-1}"],
                        curr_col,
                        is_delete=True,
                    )
            m.append(curr_col)

        # TERMINATION
        final_scores = [
            m[-1][s] + math.log(jump)
            for s in m[-1]
            if (jump := self.transitions[s].get("End", 0)) > 0
        ]
        return self._log_sum_exp(final_scores)

    def train(self, dataset):
        """Performs Viterbi training on a dataset."""
        paths = [self.decode(gene) for gene in dataset]
        self._retrain_emissions(dataset, paths)
        self._retrain_transitions(paths)

    def _retrain_emissions(self, dataset, paths):
        tally = defaultdict(lambda: defaultdict(float))
        for seq, path in zip(dataset, paths):
            char_idx = 0
            for state in path:
                if not state.startswith("D"):
                    tally[state][seq[char_idx]] += 1
                    char_idx += 1

        all_states = [f"I{k}" for k in range(self.num_matches + 1)] + [
            f"M{k}" for k in range(1, self.num_matches + 1)
        ]
        self.emissions = {
            s: self._normalize(tally[s], ALPHABET) for s in all_states
        }

    def _retrain_transitions(self, paths):
        tally = self._get_empty_transitions()
        for path in paths:
            prev = "Start"
            for state in path:
                tally[prev][state] += 1
                prev = state
            tally[prev]["End"] += 1

        self.transitions = {
            s: self._normalize(counts) for s, counts in tally.items()
        }

    def _get_empty_transitions(self):
        """Pre-initializes the strict HMM topology."""
        trans = {"Start": {"M1": 0.0, "D1": 0.0, "I0": 0.0}}
        for k in range(self.num_matches + 1):
            trans[f"I{k}"] = {f"I{k}": 0.0}
            if k < self.num_matches:
                trans[f"I{k}"].update({f"M{k+1}": 0.0, f"D{k+1}": 0.0})
            else:
                trans[f"I{k}"]["End"] = 0.0

            if k > 0:
                trans[f"M{k}"] = {f"I{k}": 0.0}
                trans[f"D{k}"] = {f"I{k}": 0.0}
                if k < self.num_matches:
                    for s in [f"M{k}", f"D{k}"]:
                        trans[s].update({f"M{k+1}": 0.0, f"D{k+1}": 0.0})
                else:
                    trans[f"M{k}"]["End"] = 0.0
                    trans[f"D{k}"]["End"] = 0.0
        return trans

    def _normalize(self, tally, alphabet=None):
        """Converts tallies to probabilities with pseudocounts."""
        pseudocount = 0.01
        keys = alphabet if alphabet is not None else tally.keys()
        total = sum(tally.values()) + (len(keys) * pseudocount)
        return {k: (tally.get(k, 0) + pseudocount) / total for k in keys}

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
        for col in zip(*msa):
            gaps = col.count("-")
            states_seq.append(MATCH if gaps < len(msa) / 2 else INSERT)
        return states_seq

    def _calculate_emissions_matrix(self, msa, states_seq):
        tally, match_cnt = defaultdict(lambda: defaultdict(float)), 0
        for i, s_type in enumerate(states_seq):
            if s_type == MATCH:
                match_cnt += 1
            key = f"{s_type}{match_cnt}"
            for seq in msa:
                if (char := seq[i]) in ALPHABET:
                    tally[key][char] += 1

        all_states = [f"I{k}" for k in range(self.num_matches + 1)] + [
            f"M{k}" for k in range(1, self.num_matches + 1)
        ]
        return {s: self._normalize(tally[s], ALPHABET) for s in all_states}

    def _calculate_transitions_matrix(self, msa, states_seq):
        tally = self._get_empty_transitions()
        for seq in msa:
            prev, match_cnt = "Start", 0
            for i, char in enumerate(seq):
                if states_seq[i] == MATCH:
                    match_cnt += 1

                if char != "-":
                    curr = f"{states_seq[i]}{match_cnt}"
                elif states_seq[i] == MATCH:
                    curr = f"D{match_cnt}"
                else:
                    continue

                tally[prev][curr] += 1
                prev = curr
            tally[prev]["End"] += 1
        return {s: self._normalize(counts) for s, counts in tally.items()}
