from src.sequence_utils import normalize_reference, normalize_sequence, reverse_complement, find_all_exact


def test_fasta_normalization_and_validation():
    seq, error = normalize_reference(">chr1\n acg t\nNN\n")
    assert error is None
    assert seq == "ACGTNN"
    _, error = normalize_reference("ACGX")
    assert "invalid base" in error


def test_reverse_complement_and_multiple_matches():
    assert reverse_complement("ACGTN") == "NACGT"
    assert find_all_exact("AAACGTTACGTT", "ACGTT") == [2, 7]
