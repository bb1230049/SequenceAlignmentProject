import random
import sys

def generate_dna_sequence(length):
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choice(bases) for _ in range(length))


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python generate_sequences.py <length> <output_prefix>")
        sys.exit()

    length = int(sys.argv[1])
    prefix = sys.argv[2]

    seq1 = generate_dna_sequence(length)
    seq2 = generate_dna_sequence(length)

    with open(prefix + "_seq1.fasta", "w") as f:
        f.write(">Sequence1\n")
        f.write(seq1 + "\n")

    with open(prefix + "_seq2.fasta", "w") as f:
        f.write(">Sequence2\n")
        f.write(seq2 + "\n")

    print("Sequences generated successfully!")