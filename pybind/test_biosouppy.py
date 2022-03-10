from biosouppy import NucleicAcid
from copy import deepcopy


class TestNucleicAcid:
  def test_code(self):
    seq = NucleicAcid(
        "test", "AaAaCcCcGgGgTtTt------ACGTRYKMSWBDHVN-nvhdbwsmkyrtgca------tTtTgGgGcCcCaAaA")

    assert 0 == seq.code(16)
    assert 1 == seq.code(23)
    assert 2 == seq.code(35)
    assert 3 == seq.code(59)

    assert "AAAACCCCGGGGTTTTAAAAAAACGTATGCCACATGAAAGTACACCGTATGCAAAAAAATTTTGGGGCCCCAAAA" \
        == seq.inflate_data()

    assert "TATGCCACATGAAAGTACACCGTAT" == seq.inflate_data(25, 25)
    assert "TGAAAGT" == seq.inflate_data(34, 7)
    assert "C" == seq.inflate_data(29, 1)
    assert "G" == seq.inflate_data(64, 1)
    assert "CCAAAA" == seq.inflate_data(69)

  def test_quality(self):
    seq = NucleicAcid("test",
                      "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",
                      "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")

    assert 31 == seq.score(42)
    assert 78 == seq.score(84)

    assert "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@oooooooooooooooooooooooooooooo" == seq.inflate_quality()

    assert "@@@@@@@@@@@@@@@@" == seq.inflate_quality(0, 16)
    assert "@@oo" == seq.inflate_quality(62, 4)
    assert "@" == seq.inflate_quality(63, 1)
    assert "o" == seq.inflate_quality(64, 1)
    assert "o" == seq.inflate_quality(93)
    assert "ooo" == seq.inflate_quality(91, 20)
    assert "" == seq.inflate_quality(95)

  def test_reverse_and_complement(self):
    seq = NucleicAcid(
        "test",
        "ACGTACTGAGCTAGTCATCGATGCCAGTCATGCGATCGTACTAGCTGAGACTGATCGCATGCTAGTACGTCA",
        "0123456789012345678901234567890123456789012345678901234567890123ZZZZZZZZ")

    cseq = deepcopy(seq)
    cseq.reverse_and_complement()

    assert "TGACGTACTAGCATGCGATCAGTCTCAGCTAGTACGATCGCATGACTGGCATCGATGACTAGCTCAGTACGT" == cseq.inflate_data()
    assert "ZZZZZZZZ4444444444444444444444444444444444444444444444444444444444444444" == cseq.inflate_quality()

    cseq.reverse_and_complement()
    assert cseq.inflate_data() == seq.inflate_data()
    assert cseq.inflate_quality() == seq.inflate_quality()
