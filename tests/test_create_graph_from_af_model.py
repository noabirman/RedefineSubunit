# tests/test_create_graph_from_af_model.py
import json
import os
import itertools
import numpy as np
import pytest
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import create_graph_from_af_model as cg
from create_graph_from_af_model import SubunitInfo

class DummySeqIO:
    """A dummy SeqIO.parse replacement that yields SeqRecords."""
    @staticmethod
    def parse(path, fmt):
        # ignore path/fmt, always yield two fragments
        yield SeqRecord(Seq("AAA"), id="seg1")
        yield SeqRecord(Seq("BBBB"), id="seg2")

@pytest.fixture(autouse=True)
def patch_seqio(monkeypatch):
    """Monkeypatch SeqIO.parse for extract_sequence_with_seqio tests."""
    monkeypatch.setattr("create_graph_from_af_model.SeqIO", DummySeqIO)
    yield

def test_extract_sequence_with_seqio_concatenates():
    # Regardless of file fmt/path, should get "AAA" + "BBBB"
    seq = cg.extract_sequence_with_seqio("dummy.cif", af_version=3)
    assert seq == "AAABBBB"

@pytest.mark.parametrize("plddt,gap,expected", [
    ([], 3, []),
    ([10,50,50,50, 10,50,50, 10,50,50,50,50,50], 2,
     # indices >40 are 1,2,3 then 5,6 then 8,9,10,11,12
     [ (1,3), (5,6), (8,12) ]),
    ([45,46,47,48,49,50], 1, [(0,5)]),
])
def test_find_high_confidence_regions(plddt, gap, expected):
    arr = np.array(plddt)
    regions = cg.find_high_confidence_regions(arr, confidence_threshold=40, gap_threshold=gap)
    assert regions == expected

def test_extract_subunit_info_single_chain():
    # A single region [start,end] on one chain
    idxs = [(2,5)]
    # token_chain_ids length at least end+1
    token_chain_ids = list("AAAABBBBB")
    full_seq = "ABCDEFGHJ"  # length 10
    infos = cg.extract_subunit_info(idxs, token_chain_ids, full_seq)
    # chain 'A' occurs positions 2-3, chain 'B' positions 4-5
    # Expect two SubunitInfos: A_1 (pos2-3) and B_1 (4-5)
    names = {si.name for si in infos}
    assert names == {"A_1", "B_1"}
    # Check sequences
    seqs = {si.name: si.sequence for si in infos}
    assert seqs["A_1"] == full_seq[2:3+1]  # "C","D"
    assert seqs["B_1"] == full_seq[4:5+1]  # "E","F"

def test_extract_subunit_info_multiple_occurrences():
    # Two regions on same chain 'X'
    idxs = [(0,2), (4,6)]
    token_chain_ids = list("XXXYYYZZ")
    full_seq = "ABCDEFGHIJ"  # len 10
    infos = cg.extract_subunit_info(idxs, token_chain_ids, full_seq)
    # Both on chain X → names X_1 and X_2
    names = [si.name for si in infos]
    assert set(names) == {"X_1", "Y_1", "Z_1"}  # Y,Z from second region
    # ordering: first region X_1, second Y_1,Z_1
    assert any(si.name == "X_1" and si.start == 0 and si.end == 2 for si in infos)

def test_find_edges_threshold():
    # build three subunits with non-overlapping indices
    s1 = SubunitInfo("A", ["A"], 0, 1, "AA")
    s2 = SubunitInfo("B", ["B"], 2, 3, "BB")
    subs = [s1, s2]
    # PAE matrix: low values => edge, high => no edge
    pae = np.array([
        [ 0, 10, 100, 100],
        [10,  0, 100, 100],
        [100,100,  0,  10],
        [100,100, 10,   0],
    ])
    # threshold=50 => edges between A–B since mean of block (0:2,2:4) is 100
    edges = cg.find_edges(subs, pae, threshold=50)
    # since 100 > threshold, not added; threshold=150 would add
    assert edges == []
    edges2 = cg.find_edges(subs, pae, threshold=150)
    assert edges2 == [("A", "B", pytest.approx(100.0))]

@pytest.mark.parametrize("filename, token_ids, expected", [
    ("A_B_full.json", list("AAAABBBB"), list("AAAA" + "BBBB")),
    ("X_Y_Z_conf.json", ["X","Y","Z","X"], ["X","Y","Z","X"]),  # only first two parts used
])
def test_rename_chains_from_file(filename, token_ids, expected):
    # create a dummy path
    path = os.path.join("/tmp", filename)
    out = cg.rename_chains_from_file(path, token_ids)
    assert out == expected
