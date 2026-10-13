"""Shared estimator specification for the regional and sensitivity analyses."""

from sklearn.ensemble import HistGradientBoostingClassifier


def regional_estimator(seed: int) -> HistGradientBoostingClassifier:
    return HistGradientBoostingClassifier(
        max_iter=300,
        max_depth=3,
        learning_rate=0.04,
        l2_regularization=1.0,
        random_state=seed,
    )
