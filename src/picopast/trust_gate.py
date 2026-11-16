"""RQ1 trust gate: one verdict per group, from a rule fixed before results are seen.

Verdicts (PRD section 3): ``usable_label``, ``relative_index_only``,
``excluded``, or ``insufficient_matchups`` (then the PRD fallback applies:
compare spatial and seasonal patterns with published cruise maps).

The thresholds live in ``config/picopast.yaml`` under ``trust_gate`` and are
recorded in ``docs/trust_gate_preregistration.md``.
"""

from __future__ import annotations

from picopast import metrics
from picopast.config import Config

VERDICTS = ("usable_label", "relative_index_only", "excluded", "insufficient_matchups")


def verdict(stats: dict[str, dict[str, float]], rule: dict) -> tuple[str, list[str]]:
    """Apply the rule to the output of ``metrics.cruise_bootstrap``.

    Returns the verdict and the list of checks that decided it.
    """
    n = stats["n"]["value"]
    if n < rule["min_matchups"]:
        return "insufficient_matchups", [f"n = {n} < {rule['min_matchups']}"]

    u = rule["usable_label"]
    checks = {
        f"log-RMSE {stats['log_rmse']['value']:.2f} <= {u['max_log_rmse']}":
            stats["log_rmse"]["value"] <= u["max_log_rmse"],
        f"|bias| {abs(stats['log_bias']['value']):.2f} <= {u['max_abs_log_bias']}":
            abs(stats["log_bias"]["value"]) <= u["max_abs_log_bias"],
        f"Spearman CI low {stats['spearman']['lo']:.2f} >= {u['min_spearman_ci_low']}":
            stats["spearman"]["lo"] >= u["min_spearman_ci_low"],
        f"within factor 2 {stats['within_factor_2']['value']:.2f} >= {u['min_within_factor_2']}":
            stats["within_factor_2"]["value"] >= u["min_within_factor_2"],
    }
    if all(checks.values()):
        return "usable_label", list(checks)

    failed = [text for text, passed in checks.items() if not passed]
    r = rule["relative_index"]
    if stats["spearman"]["value"] >= r["min_spearman"] and stats["spearman"]["lo"] > r["min_spearman_ci_low"]:
        return "relative_index_only", ["failed: " + "; ".join(failed), "rank agreement holds"]
    return "excluded", ["failed: " + "; ".join(failed), "rank agreement does not hold"]


def evaluate(cfg: Config, est_log10, obs_log10, cruise, n_boot: int = 2000, seed: int = 0) -> dict:
    """Metrics with cruise-bootstrap intervals, plus the verdict, for one group."""
    stats = metrics.cruise_bootstrap(est_log10, obs_log10, cruise, n_boot=n_boot, seed=seed)
    decision, reasons = verdict(stats, cfg["trust_gate"])
    # With one cruise the cruise bootstrap has nothing to resample.
    provisional = stats["n_cruises"]["value"] < 2 and decision != "insufficient_matchups"
    if provisional:
        reasons.append("provisional: all matchups come from one cruise")
    return {"stats": stats, "verdict": decision, "provisional": provisional, "reasons": reasons}
