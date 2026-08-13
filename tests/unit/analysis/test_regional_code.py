import pandas as pd

from regional_code_paper.analysis.regional_code import per_protein_auc


def test_per_protein_auc_keeps_full_and_composition_estimands_separate() -> None:
    tiles = pd.DataFrame(
        {
            "accession": ["P1", "P1", "P2", "P2"],
            "label": [0, 1, 0, 0],
            "full_score": [0.1, 0.9, 0.2, 0.3],
            "composition_score": [0.8, 0.2, 0.2, 0.3],
        }
    )

    result = per_protein_auc(tiles)

    assert result.to_dict("records") == [
        {
            "accession": "P1",
            "within_protein_auroc": 1.0,
            "composition_within_protein_auroc": 0.0,
        }
    ]
