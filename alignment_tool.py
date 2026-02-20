import sys

def read_fasta(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    sequence = ""
    for line in lines:
        if not line.startswith(">"):
            sequence += line.strip()
    return sequence

def read_parameters(filename):
    params = {}
    with open(filename, 'r') as f:
        for line in f:
            key, value = line.strip().split("=")
            params[key] = int(value)
    return params

def smith_waterman_affine(seq1, seq2, match, mismatch, gap_open, gap_extend):
    n = len(seq1)
    m = len(seq2)

    score = [[0]*(m+1) for _ in range(n+1)]
    gap1 = [[0]*(m+1) for _ in range(n+1)]
    gap2 = [[0]*(m+1) for _ in range(n+1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n+1):
        for j in range(1, m+1):

            if seq1[i-1] == seq2[j-1]:
                diag = score[i-1][j-1] + match
            else:
                diag = score[i-1][j-1] + mismatch

            gap1[i][j] = max(
                score[i-1][j] + gap_open,
                gap1[i-1][j] + gap_extend
            )

            gap2[i][j] = max(
                score[i][j-1] + gap_open,
                gap2[i][j-1] + gap_extend
            )

            score[i][j] = max(0, diag, gap1[i][j], gap2[i][j])

            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    align1 = ""
    align2 = ""

    i, j = max_pos

    while i > 0 and j > 0 and score[i][j] != 0:

        if seq1[i-1] == seq2[j-1]:
            score_current = match
        else:
            score_current = mismatch

        if score[i][j] == score[i-1][j-1] + score_current:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif score[i][j] == gap1[i][j]:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2, max_score

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: python alignment_tool.py seq1.fasta seq2.fasta parameters.txt")
        sys.exit()

    seq1 = read_fasta(sys.argv[1])
    seq2 = read_fasta(sys.argv[2])

    params = read_parameters(sys.argv[3])

    match = params["match"]
    mismatch = params["mismatch"]
    gap_open = params["gap_open"]
    gap_extend = params["gap_extend"]

    align1, align2, score = smith_waterman_affine(
        seq1, seq2, match, mismatch, gap_open, gap_extend
    )

    print("\nBest Local Alignment:\n")
    print(align1)
    print(align2)
    print("\nAlignment Score:", score)

# Write output to file
with open("alignment_output.txt", "w") as f:
    f.write("Best Local Alignment:\n\n")
    f.write(align1 + "\n")
    f.write(align2 + "\n\n")
    f.write("Alignment Score: " + str(score) + "\n")

print("\nOutput saved to alignment_output.txt")

    