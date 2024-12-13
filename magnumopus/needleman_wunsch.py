#!/usr/bin/env python3

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
	# Construct an empty matrix of 0's.
	matrix = []
	column_number = len(seq_a)+2
	row_number = len(seq_b)+2
	for _ in range(row_number):
		matrix.append([0]*column_number)

	# Fill in the sequences in the first row and column.
	a = 2
	b = 2
	for base_a in seq_a:
		matrix[0][a] = base_a
		a += 1
		if a == column_number:
			break
	for base_b in seq_b:
		matrix[b][0] = base_b
		b += 1
		if b == row_number:
			break

	# Fill top row and left column as gaps.
	for c in range(2, column_number):
		matrix[1][c] = (c-1)*gap
	for d in range(2, row_number):
		matrix[d][1] = (d-1)*gap

	# Fill the matrix.
	for e in range(2, row_number):
		for f in range(2, column_number):
			if matrix[0][f] == matrix[e][0]:
				match_mis_score = matrix[e-1][f-1] + match
			else:
				match_mis_score = matrix[e-1][f-1] + mismatch
			gap_score_a = matrix[e-1][f] + gap
			gap_score_b = matrix[e][f-1] + gap
			matrix[e][f] = max(match_mis_score, gap_score_a, gap_score_b)

	# Trace back the final alignment.
	seq_a_aligned = []
	seq_b_aligned = []
	g = row_number - 1
	h = column_number - 1
	score = matrix[g][h]
	while g > 1 or h > 1:
		if g > 1 and matrix[g][h] == matrix[g-1][h] + gap:
			seq_a_aligned.insert(0, "-")
			seq_b_aligned.insert(0, seq_b[g-2])
			g -= 1
		elif h > 1 and matrix[g][h] == matrix[g][h-1] + gap:
			seq_a_aligned.insert(0, seq_a[h-2])
			seq_b_aligned.insert(0, "-")
			h -= 1
		else:
			seq_a_aligned.insert(0, seq_a[h-2])
			seq_b_aligned.insert(0, seq_b[g-2])
			g -= 1
			h -= 1
	aligned_1 = "".join(seq_a_aligned)
	aligned_2 = "".join(seq_b_aligned)
	return (aligned_1, aligned_2), score