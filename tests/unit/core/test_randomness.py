from regional_code_paper.core.randomness import stable_seed


def test_stable_seed_is_deterministic_and_token_specific() -> None:
    assert stable_seed(42, "analysis", "fold", 1) == stable_seed(42, "analysis", "fold", 1)
    assert stable_seed(42, "analysis", "fold", 1) != stable_seed(42, "analysis", "fold", 2)
