from regional_code_paper.analysis.detection_aware_clustering import (
    mapped_spans,
    tryptic_span,
)


def test_tryptic_span_honors_proline_exception() -> None:
    sequence = "AAKPAARKTT"
    assert tryptic_span(sequence, 4) == (1, 7)
    assert tryptic_span(sequence, 9) == (9, 10)


def test_peptide_mapping_requires_the_site_inside_the_occurrence() -> None:
    sequence = "ASTAAST"
    assert mapped_spans(sequence, "AST", 2) == {(1, 3)}
    assert mapped_spans(sequence, "AST", 6) == {(5, 7)}
    assert mapped_spans(sequence, "A[+80]ST", 2) == set()
