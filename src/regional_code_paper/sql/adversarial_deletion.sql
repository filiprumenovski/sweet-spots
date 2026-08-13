-- Figure 1F: strongest-effect-first deletion sensitivity.
--
-- The estimand is unchanged from the primary analysis: the mean observed
-- close-pair fraction divided by the mean null fraction. `effect` determines
-- deletion order only; it is never averaged into the plotted fold.

CREATE OR REPLACE TEMP TABLE base AS
SELECT
    accession,
    observed_close_pair_fraction AS observed_fraction,
    null_close_pair_fraction AS null_fraction,
    effect,
    row_number() OVER (ORDER BY effect DESC, accession) AS adversarial_rank,
    count(*) OVER () AS n_total
FROM read_csv_auto(getenv('PAPER_CLUSTERING_INPUT'), header = true)
WHERE universe = 'all_eligible_residues'
  AND ptm = 'oglcnac'
  AND radius = 10;

CREATE OR REPLACE TEMP TABLE percentage_grid AS
SELECT
    removed_pct,
    floor((SELECT max(n_total) FROM base) * removed_pct / 100.0)::BIGINT AS removed_n
FROM range(0, 81) AS grid(removed_pct);

CREATE OR REPLACE TEMP TABLE adversarial_curve AS
SELECT
    grid.removed_pct,
    grid.removed_n,
    count(*) AS remaining_n,
    avg(base.observed_fraction) / avg(base.null_fraction) AS group_fold,
    avg(base.effect) AS mean_effect,
    sum((base.effect > 0)::INTEGER) AS positive_remaining
FROM percentage_grid AS grid
JOIN base ON base.adversarial_rank > grid.removed_n
GROUP BY ALL
ORDER BY grid.removed_pct;

CREATE OR REPLACE TEMP TABLE random_ranked AS
SELECT
    replicate.replicate,
    base.*,
    row_number() OVER (
        PARTITION BY replicate.replicate
        ORDER BY hash(base.accession || ':' || replicate.replicate::VARCHAR), base.accession
    ) AS random_rank
FROM range(1, 501) AS replicate(replicate)
CROSS JOIN base;

CREATE OR REPLACE TEMP TABLE random_curves AS
SELECT
    ranked.replicate,
    grid.removed_pct,
    grid.removed_n,
    count(*) AS remaining_n,
    avg(ranked.observed_fraction) / avg(ranked.null_fraction) AS group_fold
FROM percentage_grid AS grid
JOIN random_ranked AS ranked ON ranked.random_rank > grid.removed_n
GROUP BY ALL;

CREATE OR REPLACE TEMP TABLE random_envelope AS
SELECT
    removed_pct,
    removed_n,
    min(remaining_n) AS remaining_n,
    quantile_cont(group_fold, 0.025) AS fold_q025,
    quantile_cont(group_fold, 0.25) AS fold_q25,
    quantile_cont(group_fold, 0.5) AS fold_median,
    quantile_cont(group_fold, 0.75) AS fold_q75,
    quantile_cont(group_fold, 0.975) AS fold_q975
FROM random_curves
GROUP BY ALL
ORDER BY removed_pct;

CREATE OR REPLACE TEMP TABLE single_step_curve AS
WITH steps AS (
    SELECT range AS removed_n
    FROM range(0, (SELECT max(n_total) FROM base))
)
SELECT
    steps.removed_n,
    100.0 * steps.removed_n / max(base.n_total) AS removed_pct,
    count(*) AS remaining_n,
    avg(base.observed_fraction) / avg(base.null_fraction) AS group_fold
FROM steps
JOIN base ON base.adversarial_rank > steps.removed_n
GROUP BY steps.removed_n
ORDER BY steps.removed_n;

COPY adversarial_curve
TO (getenv('PAPER_OUTPUT_ROOT') || '/adversarial_curve.csv')
(HEADER, DELIMITER ',');

COPY random_envelope
TO (getenv('PAPER_OUTPUT_ROOT') || '/random_envelope.csv')
(HEADER, DELIMITER ',');

COPY (
    WITH crossing AS (
        SELECT min(removed_n) FILTER (WHERE group_fold < 1) AS first_below_null_n
        FROM single_step_curve
    )
    SELECT
        count(*) AS n_proteins,
        sum((effect > 0)::INTEGER) AS n_positive_effect,
        avg(observed_fraction) / avg(null_fraction) AS full_group_fold,
        (SELECT first_below_null_n FROM crossing) AS first_below_null_n,
        100.0 * (SELECT first_below_null_n FROM crossing) / count(*) AS first_below_null_pct,
        (SELECT group_fold FROM single_step_curve WHERE removed_n = 1) AS fold_after_top_one,
        (SELECT group_fold FROM adversarial_curve WHERE removed_pct = 10) AS fold_after_top_10pct,
        (SELECT group_fold FROM adversarial_curve WHERE removed_pct = 25) AS fold_after_top_25pct,
        (SELECT group_fold FROM adversarial_curve WHERE removed_pct = 50) AS fold_after_top_50pct,
        (SELECT group_fold FROM adversarial_curve WHERE removed_pct = 70) AS fold_after_top_70pct
    FROM base
)
TO (getenv('PAPER_OUTPUT_ROOT') || '/summary.csv')
(HEADER, DELIMITER ',');
