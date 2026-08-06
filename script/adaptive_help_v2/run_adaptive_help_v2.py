#!/usr/bin/env python3
"""
PACS adaptive HELP normalization v2

目标：
1. 不使用 HC/T2DM 标签计算任何 normalization factor。
2. HELP 分别用于仪器、板级 detection、样本级异常判断。
3. pooled QC endogenous proteins 用于补充板级前处理校正。
4. 比较多个候选流程，以技术偏差移除为主、组内聚集/组间分离为验证。
"""

from __future__ import annotations

import json
import math
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, silhouette_score
from sklearn.model_selection import LeaveOneGroupOut, RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=RuntimeWarning)
SEED = 20260806
np.random.seed(SEED)

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = ROOT / "Results" / "adaptive_help_v2"
for sub in ["audit", "factors", "metrics", "matrices", "plots"]:
    (OUT / sub).mkdir(parents=True, exist_ok=True)

OUTLIERS = {
    "WJ-1-5","WJ-1-7","WJ-1-19","WJ-1-21","WJ-1-44","WJ-1-67",
    "WJ-1-90","WJ-1-101","WJ-1-102","WJ-1-103","WJ-1-104",
    "WJ-1-105","WJ-1-112","WJ-1-146","WJ-1-148","WJ-1-157",
    "WJ-1-171","WJ-1-176","WJ-1-177","WJ-1-178","WJ-1-187",
    "WJ-1-203","WJ-1-208","WJ-1-209","WJ-1-212","WJ-1-213",
    "WJ-1-43","WJ-1-59","WJ-1-83",
}


def robust_mad(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return np.nan
    return float(np.median(np.abs(x - np.median(x))))


def robust_z(x: pd.Series) -> pd.Series:
    med = float(np.nanmedian(x))
    scale = 1.4826 * robust_mad(x.to_numpy())
    if not np.isfinite(scale) or scale < 1e-8:
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - med) / scale


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if ok.sum() == 0:
        return np.nan
    values, weights = values[ok], weights[ok]
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cutoff = 0.5 * weights.sum()
    return float(values[np.searchsorted(np.cumsum(weights), cutoff)])


def standardize_metadata(path: Path, study: bool) -> pd.DataFrame:
    x = pd.read_excel(path)
    x.columns = [str(c).strip() for c in x.columns]
    key = {str(c).lower().replace(" ", "_"): c for c in x.columns}

    def pick(*names, required=True):
        for name in names:
            if name in key:
                return key[name]
        for normalized, original in key.items():
            if any(name in normalized for name in names):
                return original
        if required:
            raise KeyError(f"Missing metadata column: {names}")
        return None

    sid = pick("sample_id", "sampleid")
    sname = pick("sample_name", "samplename", required=False)
    group1 = pick("group1", "group_1", required=False)
    group2 = pick("group2", "group_2", required=False)
    plate = pick("plate", "plate_id", required=False)
    run = pick("run", "run_id", "instrument_run", required=False)

    out = pd.DataFrame({
        "sample_id": x[sid].astype(str).str.strip(),
        "sample_name": x[sname].astype(str).str.strip() if sname else "",
        "group1": x[group1].astype(str).str.strip() if group1 else "",
        "group2": x[group2].astype(str).str.strip() if group2 else "",
        "Plate": x[plate].astype(str).str.strip() if plate else "",
        "Run": x[run].astype(str).str.strip() if run else "",
    })
    out = out.replace({"nan": "", "None": ""})
    out["group2"] = out["group2"].where(out["group2"].ne(""), out["group1"])

    if study:
        out["sample_role"] = "study"
    else:
        all_text = x.fillna("").astype(str).agg(" | ".join, axis=1).str.lower()
        role = np.select(
            [
                all_text.str.contains("blank", regex=False),
                all_text.str.contains("neat", regex=False) | all_text.str.contains("qc2", regex=False),
                all_text.str.contains("internal|calibrator|qc3", regex=True),
                all_text.str.contains("qc1|plate.?qc", regex=True),
            ],
            ["blank", "neat", "internal_qc", "plate_qc"],
            default="other_control",
        )
        out["sample_role"] = role
    return out


def read_help(path: Path, known_ids: set[str]) -> pd.DataFrame:
    x = pd.read_excel(path)
    first = x.columns[0]
    col_overlap = sum(str(c) in known_ids for c in x.columns[1:])
    row_overlap = x[first].astype(str).isin(known_ids).sum()

    if col_overlap >= row_overlap:
        m = x.set_index(first)
    else:
        m = x.set_index(first).T
    m.index = m.index.astype(str)
    m.columns = m.columns.astype(str)
    m = m.apply(pd.to_numeric, errors="coerce")
    m[m <= 0] = np.nan
    raw_med = np.nanmedian(m.to_numpy())
    if raw_med > 100:
        m = np.log2(m)
    return m


def factor_from_help(help_mat: pd.DataFrame, ids: list[str], features: list[str],
                     reference: pd.Series, weights: pd.Series) -> pd.Series:
    ids = [x for x in ids if x in help_mat.columns]
    f = [x for x in features if x in help_mat.index]
    dev = help_mat.loc[f, ids].sub(reference.loc[f], axis=0)
    out = {}
    for sid in ids:
        out[sid] = weighted_median(dev[sid].to_numpy(), weights.loc[f].to_numpy())
    return pd.Series(out, dtype=float)


def grouped_shrink(factors: pd.Series, meta: pd.DataFrame, group_col: str,
                   min_effect: float = 0.03) -> tuple[pd.DataFrame, float]:
    z = meta[["sample_id", group_col]].copy()
    z["factor"] = z["sample_id"].map(factors)
    z = z.dropna(subset=["factor"])
    grouped = z.groupby(group_col)["factor"].agg(["median", "var", "count"]).reset_index()
    grouped["var"] = grouped["var"].fillna(0.0)
    between = float(np.nanvar(grouped["median"], ddof=1)) if len(grouped) > 1 else 0.0
    noise = float(np.nanmedian(grouped["var"] / grouped["count"].clip(lower=1)))
    reliability = 0.0 if between <= 1e-12 else float(np.clip((between - noise) / between, 0, 1))
    if grouped["median"].abs().max() < min_effect:
        reliability = 0.0
    grouped["centered_factor"] = grouped["median"] - grouped["median"].median()
    grouped["reliability"] = reliability
    grouped["applied_factor"] = grouped["centered_factor"] * reliability
    return grouped, reliability


def map_group_factor(meta: pd.DataFrame, factor_table: pd.DataFrame, group_col: str) -> pd.Series:
    lookup = factor_table.set_index(group_col)["applied_factor"]
    return meta.set_index("sample_id")[group_col].map(lookup).fillna(0.0)


def apply_scalar(mat: pd.DataFrame, scalar: pd.Series) -> pd.DataFrame:
    scalar = scalar.reindex(mat.columns).fillna(0.0)
    return mat.sub(scalar, axis=1)


def estimate_protein_plate_offsets(qc_mat: pd.DataFrame, qc_meta: pd.DataFrame
                                   ) -> tuple[pd.DataFrame, pd.Series]:
    plates = sorted(qc_meta["Plate"].dropna().astype(str).unique())
    global_ref = qc_mat.median(axis=1, skipna=True)
    offsets = pd.DataFrame(index=qc_mat.index, columns=plates, dtype=float)
    within_vars = []
    counts = []
    for plate in plates:
        ids = qc_meta.loc[qc_meta["Plate"].astype(str).eq(plate), "sample_id"]
        ids = [x for x in ids if x in qc_mat.columns]
        sub = qc_mat.loc[:, ids]
        offsets[plate] = sub.median(axis=1, skipna=True) - global_ref
        within_vars.append(sub.var(axis=1, skipna=True))
        counts.append(sub.notna().sum(axis=1).clip(lower=1))
    within_var = pd.concat(within_vars, axis=1).median(axis=1, skipna=True)
    n_eff = pd.concat(counts, axis=1).median(axis=1, skipna=True)
    between_var = offsets.var(axis=1, skipna=True)
    noise = within_var / n_eff
    reliability = ((between_var - noise) / between_var.replace(0, np.nan)).clip(0, 1).fillna(0)
    detect = qc_mat.notna().mean(axis=1)
    reliability[(detect < 0.80) | (offsets.abs().max(axis=1) < 0.05)] = 0.0
    applied = offsets.mul(reliability, axis=0)
    return applied, reliability


def apply_protein_plate(mat: pd.DataFrame, meta: pd.DataFrame,
                        offsets: pd.DataFrame) -> pd.DataFrame:
    out = mat.copy()
    plate_lookup = meta.set_index("sample_id")["Plate"].astype(str)
    for sid in out.columns:
        plate = plate_lookup.get(sid, "")
        if plate in offsets.columns:
            out[sid] = out[sid] - offsets[plate]
    return out


def sample_help_decisions(help_mat: pd.DataFrame, study_ids: list[str],
                          panel: list[str], reference: pd.Series, weights: pd.Series,
                          run_scalar: pd.Series, plate_scalar: pd.Series,
                          protein_mat: pd.DataFrame) -> pd.DataFrame:
    f = [x for x in panel if x in help_mat.index]
    dev = help_mat.loc[f, study_ids].sub(reference.loc[f], axis=0)
    technical = run_scalar.reindex(study_ids).fillna(0) + plate_scalar.reindex(study_ids).fillna(0)
    dev = dev.sub(technical, axis=1)

    rows = []
    cohort_protein_median = protein_mat.median(axis=0, skipna=True)
    center_protein = cohort_protein_median.median()
    for sid in study_ids:
        vals = dev[sid].to_numpy()
        factor = weighted_median(vals, weights.loc[f].to_numpy())
        ok = np.isfinite(vals)
        coherence = np.nan
        if ok.sum() > 0 and np.isfinite(factor) and abs(factor) > 1e-8:
            coherence = float(np.mean(np.sign(vals[ok]) == np.sign(factor)))
        protein_shift = float(cohort_protein_median.get(sid, np.nan) - center_protein)
        concordant = bool(
            np.isfinite(factor) and np.isfinite(protein_shift)
            and abs(factor) >= 0.03 and abs(protein_shift) >= 0.03
            and np.sign(factor) == np.sign(protein_shift)
        )
        rows.append({
            "sample_id": sid,
            "sample_help_factor": factor,
            "help_coherence": coherence,
            "protein_global_shift": protein_shift,
            "help_protein_concordant": concordant,
            "n_help": int(ok.sum()),
        })
    out = pd.DataFrame(rows)
    out["help_z"] = robust_z(out["sample_help_factor"])
    out["alpha"] = 0.0
    mid = (
        out["help_z"].abs().between(1.5, 2.5, inclusive="left")
        & out["help_coherence"].ge(0.65)
        & out["help_protein_concordant"]
        & out["n_help"].ge(15)
    )
    high = (
        out["help_z"].abs().ge(2.5)
        & out["help_coherence"].ge(0.70)
        & out["help_protein_concordant"]
        & out["n_help"].ge(15)
    )
    out.loc[mid, "alpha"] = 0.25
    out.loc[high, "alpha"] = 0.50
    out["applied_sample_factor"] = out["alpha"] * out["sample_help_factor"]
    out["decision"] = "no_correction"
    out.loc[out["alpha"].eq(0.25), "decision"] = "partial_25pct"
    out.loc[out["alpha"].eq(0.50), "decision"] = "partial_50pct"
    out.loc[out["help_z"].abs().ge(3.5) & out["alpha"].eq(0), "decision"] = "flag_extreme_not_concordant"
    return out


def eta_squared(scores: np.ndarray, groups: pd.Series, weights: np.ndarray) -> float:
    vals = []
    groups = groups.astype(str).to_numpy()
    for j in range(scores.shape[1]):
        y = scores[:, j]
        grand = np.mean(y)
        ss_total = np.sum((y - grand) ** 2)
        if ss_total <= 0:
            vals.append(0.0)
            continue
        ss_between = 0.0
        for g in np.unique(groups):
            yg = y[groups == g]
            ss_between += len(yg) * (np.mean(yg) - grand) ** 2
        vals.append(ss_between / ss_total)
    w = weights[:len(vals)]
    return float(np.average(vals, weights=w)) if np.sum(w) > 0 else float(np.mean(vals))


def prepare_pca(mat: pd.DataFrame, feature_ids: list[str], n_components: int = 10):
    x = mat.loc[feature_ids].T
    x = SimpleImputer(strategy="median").fit_transform(x)
    x = StandardScaler().fit_transform(x)
    n_components = min(n_components, x.shape[0] - 1, x.shape[1])
    pca = PCA(n_components=n_components, random_state=SEED)
    scores = pca.fit_transform(x)
    return scores, pca.explained_variance_ratio_


def group_metrics(scores: np.ndarray, labels: pd.Series) -> tuple[float, float, float]:
    labels = labels.astype(str).to_numpy()
    centroids = {g: scores[labels == g].mean(axis=0) for g in np.unique(labels)}
    within = np.mean([
        np.linalg.norm(scores[i] - centroids[labels[i]])
        for i in range(len(labels))
    ])
    gs = list(centroids)
    if len(gs) == 2:
        between = float(np.linalg.norm(centroids[gs[0]] - centroids[gs[1]]))
    else:
        d = []
        for i in range(len(gs)):
            for j in range(i + 1, len(gs)):
                d.append(np.linalg.norm(centroids[gs[i]] - centroids[gs[j]]))
        between = float(np.mean(d))
    ratio = between / within if within > 0 else np.nan
    sil = float(silhouette_score(scores, labels)) if len(np.unique(labels)) > 1 else np.nan
    return float(within), ratio, sil


def cv_auc(mat: pd.DataFrame, meta: pd.DataFrame, n_features: int = 20
           ) -> tuple[float, float, float]:
    y = meta["group2"].astype(str)
    valid_groups = sorted(y.value_counts().index[:2])
    keep = y.isin(valid_groups)
    y = (y[keep] == valid_groups[-1]).astype(int).to_numpy()
    ids = meta.loc[keep, "sample_id"].tolist()
    x = mat.loc[:, ids].T.to_numpy()
    plates = meta.loc[keep, "Plate"].astype(str).to_numpy()

    composite = meta.loc[keep, "group2"].astype(str) + "_P" + meta.loc[keep, "Plate"].astype(str)
    min_stratum = int(composite.value_counts().min())
    n_splits = max(2, min(3, min_stratum))
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=5, random_state=SEED)
    aucs = []
    for train, test in cv.split(x, composite):
        imp = SimpleImputer(strategy="median")
        xtr = imp.fit_transform(x[train])
        xte = imp.transform(x[test])
        finite_var = np.nanvar(xtr, axis=0) > 1e-10
        xtr, xte = xtr[:, finite_var], xte[:, finite_var]
        if xtr.shape[1] == 0:
            continue
        f, _ = f_classif(xtr, y[train])
        f = np.nan_to_num(f, nan=-np.inf)
        top = np.argsort(f)[::-1][:min(n_features, xtr.shape[1])]
        scaler = StandardScaler()
        xtr = scaler.fit_transform(xtr[:, top])
        xte = scaler.transform(xte[:, top])
        model = LogisticRegression(
            penalty="l2", C=0.2, class_weight="balanced",
            solver="liblinear", max_iter=2000, random_state=SEED
        )
        model.fit(xtr, y[train])
        aucs.append(roc_auc_score(y[test], model.predict_proba(xte)[:, 1]))

    logo_aucs = []
    for train, test in LeaveOneGroupOut().split(x, y, plates):
        if len(np.unique(y[test])) < 2 or len(np.unique(y[train])) < 2:
            continue
        imp = SimpleImputer(strategy="median")
        xtr = imp.fit_transform(x[train])
        xte = imp.transform(x[test])
        finite_var = np.nanvar(xtr, axis=0) > 1e-10
        xtr, xte = xtr[:, finite_var], xte[:, finite_var]
        f, _ = f_classif(xtr, y[train])
        f = np.nan_to_num(f, nan=-np.inf)
        top = np.argsort(f)[::-1][:min(n_features, xtr.shape[1])]
        scaler = StandardScaler()
        xtr = scaler.fit_transform(xtr[:, top])
        xte = scaler.transform(xte[:, top])
        model = LogisticRegression(
            penalty="l2", C=0.2, class_weight="balanced",
            solver="liblinear", max_iter=2000, random_state=SEED
        )
        model.fit(xtr, y[train])
        logo_aucs.append(roc_auc_score(y[test], model.predict_proba(xte)[:, 1]))

    return float(np.mean(aucs)), float(np.std(aucs)), float(np.mean(logo_aucs))


def qc_median_sd(qc_mat: pd.DataFrame) -> float:
    return float(qc_mat.std(axis=1, skipna=True).median())


def save_pca_plot(mat: pd.DataFrame, feature_ids: list[str], meta: pd.DataFrame,
                  path: Path, title: str):
    scores, var = prepare_pca(mat, feature_ids, n_components=2)
    fig, ax = plt.subplots(figsize=(7, 5))
    for group in sorted(meta["group2"].astype(str).unique()):
        idx = meta["group2"].astype(str).eq(group).to_numpy()
        ax.scatter(scores[idx, 0], scores[idx, 1], s=28, alpha=0.75, label=group)
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}%)")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    study_meta = standardize_metadata(DATA / "metadata" / "all_sample_metadata.xlsx", True)
    qc_meta = standardize_metadata(DATA / "metadata" / "all_QC_metadata.xlsx", False)
    study_meta = study_meta[~study_meta["sample_id"].isin(OUTLIERS)].copy()
    controls = qc_meta[qc_meta["sample_role"].isin(["internal_qc", "plate_qc"])].copy()

    pg = pd.read_csv(DATA / "pg_matrix" / "pg_matrix.tsv", sep="\t")
    annotation = pg.iloc[:, :5].copy()
    protein = pg.iloc[:, 5:].apply(pd.to_numeric, errors="coerce")
    protein.index = pg.iloc[:, 0].astype(str)
    protein[protein <= 0] = np.nan
    protein = np.log2(protein)

    known = set(study_meta["sample_id"]) | set(qc_meta["sample_id"])
    help_mat = read_help(DATA / "HELP" / "ALL-HELP.xlsx", known)

    study_ids = [
        x for x in study_meta["sample_id"]
        if x in protein.columns and x in help_mat.columns
    ]
    internal_ids = [
        x for x in controls.loc[controls["sample_role"].eq("internal_qc"), "sample_id"]
        if x in protein.columns and x in help_mat.columns
    ]
    plate_qc_ids = [
        x for x in controls.loc[controls["sample_role"].eq("plate_qc"), "sample_id"]
        if x in protein.columns and x in help_mat.columns
    ]
    if len(internal_ids) < 3 or len(plate_qc_ids) < 3:
        raise RuntimeError("Internal-QC or pooled-QC sample count is insufficient.")

    study_meta = study_meta.set_index("sample_id").loc[study_ids].reset_index()
    internal_meta = controls.set_index("sample_id").loc[internal_ids].reset_index()
    plate_qc_meta = controls.set_index("sample_id").loc[plate_qc_ids].reset_index()
    all_use_meta = pd.concat([study_meta, internal_meta, plate_qc_meta], ignore_index=True)

    study_raw = protein.loc[:, study_ids]
    internal_raw = protein.loc[:, internal_ids]
    plate_qc_raw = protein.loc[:, plate_qc_ids]

    audit = pd.DataFrame(index=help_mat.index)
    audit["detect_study"] = help_mat.loc[:, study_ids].notna().mean(axis=1)
    audit["detect_internal"] = help_mat.loc[:, internal_ids].notna().mean(axis=1)
    audit["detect_plate_qc"] = help_mat.loc[:, plate_qc_ids].notna().mean(axis=1)
    audit["mad_study"] = help_mat.loc[:, study_ids].apply(lambda r: robust_mad(r.to_numpy()), axis=1)
    audit["mad_internal"] = help_mat.loc[:, internal_ids].apply(lambda r: robust_mad(r.to_numpy()), axis=1)
    audit["mad_plate_qc"] = help_mat.loc[:, plate_qc_ids].apply(lambda r: robust_mad(r.to_numpy()), axis=1)

    eligible = audit[
        audit["detect_study"].ge(0.98)
        & audit["detect_internal"].ge(1.0)
        & audit["detect_plate_qc"].ge(1.0)
    ].copy()
    if len(eligible) < 15:
        eligible = audit[
            audit["detect_study"].ge(0.95)
            & audit["detect_internal"].ge(0.95)
            & audit["detect_plate_qc"].ge(0.95)
        ].copy()
    for col in ["mad_internal", "mad_plate_qc", "mad_study"]:
        eligible[f"z_{col}"] = robust_z(eligible[col].fillna(eligible[col].median())).abs()
    eligible["technical_score"] = (
        eligible["z_mad_internal"]
        + eligible["z_mad_plate_qc"]
        + 0.25 * eligible["z_mad_study"]
        + 3 * (1 - eligible["detect_study"])
    )
    ranked = eligible.sort_values("technical_score").index.tolist()

    selected = []
    corr = help_mat.loc[ranked, study_ids].T.corr(min_periods=max(20, len(study_ids)//3))
    for peptide in ranked:
        if all(abs(corr.loc[peptide, prev]) < 0.98 for prev in selected if np.isfinite(corr.loc[peptide, prev])):
            selected.append(peptide)
        if len(selected) >= 25:
            break
    for peptide in ranked:
        if peptide not in selected:
            selected.append(peptide)
        if len(selected) >= min(25, len(ranked)):
            break
    if len(selected) < 5:
        raise RuntimeError("Fewer than five reliable HELP peptides were found.")

    audit["selected"] = audit.index.isin(selected)
    audit["selection_rank"] = audit.index.map({x: i + 1 for i, x in enumerate(selected)})
    audit.reset_index(names="HELP").to_csv(OUT / "audit" / "help_panel.csv", index=False)

    controls_help_ids = internal_ids + plate_qc_ids
    ref = help_mat.loc[selected, controls_help_ids].median(axis=1, skipna=True)
    technical_mad = help_mat.loc[selected, controls_help_ids].apply(
        lambda r: robust_mad(r.to_numpy()), axis=1
    ).replace(0, np.nan)
    weights = (1.0 / technical_mad.pow(2)).replace([np.inf, -np.inf], np.nan)
    weights = weights.fillna(weights.median()).clip(upper=weights.quantile(0.95))

    internal_factor = factor_from_help(help_mat, internal_ids, selected, ref, weights)
    run_table, run_reliability = grouped_shrink(internal_factor, internal_meta, "Run", min_effect=0.03)
    run_table.to_csv(OUT / "factors" / "run_factors.csv", index=False)

    run_scalar_all = map_group_factor(all_use_meta, run_table, "Run")
    run_study = run_scalar_all.reindex(study_ids).fillna(0)
    run_internal = run_scalar_all.reindex(internal_ids).fillna(0)
    run_plate_qc = run_scalar_all.reindex(plate_qc_ids).fillna(0)

    study_run = apply_scalar(study_raw, run_study)
    internal_run = apply_scalar(internal_raw, run_internal)
    plate_qc_run = apply_scalar(plate_qc_raw, run_plate_qc)

    plate_help_factor = factor_from_help(help_mat, plate_qc_ids, selected, ref, weights) - run_plate_qc
    plate_help_table, plate_help_reliability = grouped_shrink(
        plate_help_factor, plate_qc_meta, "Plate", min_effect=0.03
    )
    plate_help_table.to_csv(OUT / "factors" / "plate_help_factors.csv", index=False)
    plate_help_scalar_all = map_group_factor(all_use_meta, plate_help_table, "Plate")
    plate_help_study = plate_help_scalar_all.reindex(study_ids).fillna(0)
    plate_help_internal = plate_help_scalar_all.reindex(internal_ids).fillna(0)
    plate_help_qc = plate_help_scalar_all.reindex(plate_qc_ids).fillna(0)

    study_help_plate = apply_scalar(study_run, plate_help_study)
    internal_help_plate = apply_scalar(internal_run, plate_help_internal)
    plate_qc_help_plate = apply_scalar(plate_qc_run, plate_help_qc)

    protein_offsets, protein_reliability = estimate_protein_plate_offsets(
        plate_qc_help_plate, plate_qc_meta
    )
    protein_reliability.rename("reliability").to_csv(
        OUT / "factors" / "protein_plate_reliability.csv"
    )
    study_qc_plate = apply_protein_plate(study_help_plate, study_meta, protein_offsets)
    internal_qc_plate = apply_protein_plate(internal_help_plate, internal_meta, protein_offsets)
    plate_qc_qc_plate = apply_protein_plate(plate_qc_help_plate, plate_qc_meta, protein_offsets)

    decisions = sample_help_decisions(
        help_mat, study_ids, selected, ref, weights,
        run_study, plate_help_study, study_qc_plate
    )
    decisions.to_csv(OUT / "factors" / "sample_help_decisions.csv", index=False)
    adaptive_scalar = decisions.set_index("sample_id")["applied_sample_factor"]
    fixed25 = decisions.set_index("sample_id")["sample_help_factor"].fillna(0) * 0.25
    fixed50 = decisions.set_index("sample_id")["sample_help_factor"].fillna(0) * 0.50

    strategies = {
        "S0_raw": (study_raw, internal_raw, plate_qc_raw),
        "S1_help_run": (study_run, internal_run, plate_qc_run),
        "S2_help_run_plate": (study_help_plate, internal_help_plate, plate_qc_help_plate),
        "S3_help_plus_qc_plate": (study_qc_plate, internal_qc_plate, plate_qc_qc_plate),
        "S4_adaptive_sample_help": (
            apply_scalar(study_qc_plate, adaptive_scalar),
            internal_qc_plate,
            plate_qc_qc_plate,
        ),
        "S5_fixed_sample_help_25pct": (
            apply_scalar(study_qc_plate, fixed25),
            internal_qc_plate,
            plate_qc_qc_plate,
        ),
        "S6_fixed_sample_help_50pct": (
            apply_scalar(study_qc_plate, fixed50),
            internal_qc_plate,
            plate_qc_qc_plate,
        ),
    }

    detect_study = study_raw.notna().mean(axis=1)
    detect_qc = plate_qc_raw.notna().mean(axis=1)
    study_mad = study_raw.apply(lambda r: robust_mad(r.to_numpy()), axis=1)
    qc_mad = plate_qc_raw.apply(lambda r: robust_mad(r.to_numpy()), axis=1)
    ratio = study_mad / qc_mad.replace(0, np.nan)
    feature_table = pd.DataFrame({
        "detect_study": detect_study,
        "detect_plate_qc": detect_qc,
        "study_mad": study_mad,
        "plate_qc_mad": qc_mad,
        "biology_technical_ratio": ratio,
    })
    feature_pool = feature_table[
        feature_table["detect_study"].ge(0.70)
        & feature_table["detect_plate_qc"].ge(0.60)
        & feature_table["study_mad"].gt(0)
    ].sort_values(["biology_technical_ratio", "study_mad"], ascending=False)
    pca_features = feature_pool.head(min(500, len(feature_pool))).index.tolist()
    feature_table.assign(selected_for_pca=feature_table.index.isin(pca_features)).reset_index(
        names="Protein.Group"
    ).to_csv(OUT / "audit" / "protein_feature_audit.csv", index=False)

    raw_var = float(np.nanmedian(study_raw.loc[pca_features].var(axis=1, skipna=True)))
    raw_qc_sd = qc_median_sd(plate_qc_raw.loc[pca_features])
    rows = []
    for name, (study_m, internal_m, qc_m) in strategies.items():
        scores, explained = prepare_pca(study_m, pca_features, n_components=10)
        plate_eta = eta_squared(scores[:, :5], study_meta["Plate"], explained[:5])
        run_eta = eta_squared(scores[:, :5], study_meta["Run"], explained[:5])
        within, between_within, sil = group_metrics(scores, study_meta["group2"])
        qcsd = qc_median_sd(qc_m.loc[pca_features])
        variance = float(np.nanmedian(study_m.loc[pca_features].var(axis=1, skipna=True)))
        variance_ratio = variance / raw_var if raw_var > 0 else np.nan
        auc_mean, auc_sd, plate_auc = cv_auc(study_m.loc[pca_features], study_meta, 20)
        rows.append({
            "strategy": name,
            "qc_median_log2_sd": qcsd,
            "qc_sd_improvement_pct": 100 * (raw_qc_sd - qcsd) / raw_qc_sd,
            "plate_eta2_pc1_5": plate_eta,
            "run_eta2_pc1_5": run_eta,
            "within_group_distance": within,
            "between_within_ratio": between_within,
            "group_silhouette": sil,
            "study_variance_ratio_vs_raw": variance_ratio,
            "cv_auc_mean": auc_mean,
            "cv_auc_sd": auc_sd,
            "leave_one_plate_out_auc": plate_auc,
        })

    metrics = pd.DataFrame(rows)
    raw_row = metrics.set_index("strategy").loc["S0_raw"]
    variance_penalty = (metrics["study_variance_ratio_vs_raw"] - 1).abs()
    metrics["label_free_score"] = (
        0.35 * (raw_row["qc_median_log2_sd"] / metrics["qc_median_log2_sd"].clip(lower=1e-6))
        + 0.25 * (raw_row["plate_eta2_pc1_5"] / metrics["plate_eta2_pc1_5"].clip(lower=1e-6))
        + 0.20 * (raw_row["run_eta2_pc1_5"] / metrics["run_eta2_pc1_5"].clip(lower=1e-6))
        + 0.20 * (1 - variance_penalty.clip(0, 1))
    )
    safety = (
        metrics["plate_eta2_pc1_5"].le(raw_row["plate_eta2_pc1_5"] * 1.05)
        & metrics["run_eta2_pc1_5"].le(raw_row["run_eta2_pc1_5"] * 1.05)
        & metrics["study_variance_ratio_vs_raw"].between(0.75, 1.25)
    )
    candidates = metrics[safety].copy()
    if candidates.empty:
        candidates = metrics[metrics["strategy"].isin(["S0_raw", "S1_help_run", "S2_help_run_plate"])].copy()
    selected_strategy = candidates.sort_values("label_free_score", ascending=False).iloc[0]["strategy"]
    metrics["selected"] = metrics["strategy"].eq(selected_strategy)
    metrics.to_csv(OUT / "metrics" / "strategy_summary.csv", index=False)

    final_study, final_internal, final_qc = strategies[selected_strategy]
    normalized = pd.concat(
        [annotation.reset_index(drop=True), final_study.reset_index(drop=True)],
        axis=1
    )
    normalized.to_csv(
        OUT / "matrices" / "final_normalized_study_matrix.tsv.gz",
        sep="\t", index=False, compression="gzip"
    )

    save_pca_plot(study_raw, pca_features, study_meta, OUT / "plots" / "pca_raw.png", "Raw")
    save_pca_plot(final_study, pca_features, study_meta, OUT / "plots" / "pca_final.png",
                  f"Final: {selected_strategy}")

    fig, ax = plt.subplots(figsize=(8, 4.8))
    plot_metrics = metrics.set_index("strategy")
    plot_metrics[["plate_eta2_pc1_5", "run_eta2_pc1_5"]].plot(kind="bar", ax=ax)
    ax.set_ylabel("Technical eta² (lower is better)")
    ax.set_xlabel("")
    ax.set_title("Residual technical variation")
    fig.tight_layout()
    fig.savefig(OUT / "plots" / "technical_bias_comparison.png", dpi=160)
    plt.close(fig)

    selected_row = metrics.set_index("strategy").loc[selected_strategy]
    summary = {
        "selected_strategy": selected_strategy,
        "n_study_samples": len(study_ids),
        "n_internal_qc": len(internal_ids),
        "n_plate_qc": len(plate_qc_ids),
        "n_help_selected": len(selected),
        "selected_help": selected,
        "run_reliability": run_reliability,
        "plate_help_reliability": plate_help_reliability,
        "n_sample_partial_25": int((decisions["alpha"] == 0.25).sum()),
        "n_sample_partial_50": int((decisions["alpha"] == 0.50).sum()),
        "metrics": {k: (None if pd.isna(v) else float(v)) for k, v in selected_row.items()
                    if k not in ["selected"]},
        "raw_metrics": {k: (None if pd.isna(v) else float(v)) for k, v in raw_row.items()},
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2))

    md = f"""# Adaptive HELP normalization v2 — 运行结果

## 最终选择

**{selected_strategy}**

本次流程没有使用 HC/T2DM 标签计算任何 normalization factor。标签只用于最后检查组内聚集和组间分离。

## 主要处理

1. 去除 README 中已确认的 29 个异常样本，并对蛋白强度做 log2。
2. 从 50 条 HELP 中，按检出率、Internal-QC 稳定性、pooled-QC 稳定性和冗余程度，自动选择 {len(selected)} 条 Core HELP。
3. 用 Internal-QC HELP 判断并校正仪器 Run 偏移。
4. 用 pooled-QC HELP 校正板级 detection 偏移，再用 pooled-QC endogenous proteins 校正板级前处理偏移。
5. 每个研究样本先判断 HELP 偏移是否真实且与整体蛋白偏移一致，再决定 0%、25% 或 50% 校正；不再对所有样本一刀切。
6. 比较 7 个候选流程，以 QC 变异、Plate/Run 残留和生物变异保留作为主选择指标；HC/T2DM 聚集和 AUC 只作验证。

## 关键结果

- Study samples: **{len(study_ids)}**
- Core HELP: **{len(selected)}**
- Run correction reliability: **{run_reliability:.3f}**
- Plate HELP correction reliability: **{plate_help_reliability:.3f}**
- 25% sample HELP correction: **{int((decisions["alpha"] == 0.25).sum())} samples**
- 50% sample HELP correction: **{int((decisions["alpha"] == 0.50).sum())} samples**
- QC median log2 SD: **{selected_row["qc_median_log2_sd"]:.4f}**
- QC SD improvement: **{selected_row["qc_sd_improvement_pct"]:.1f}%**
- Plate eta²: **{selected_row["plate_eta2_pc1_5"]:.4f}**
- Run eta²: **{selected_row["run_eta2_pc1_5"]:.4f}**
- Within-group distance: **{selected_row["within_group_distance"]:.4f}**
- Between/within ratio: **{selected_row["between_within_ratio"]:.4f}**
- Repeated CV AUC: **{selected_row["cv_auc_mean"]:.3f} ± {selected_row["cv_auc_sd"]:.3f}**
- Leave-one-plate-out AUC: **{selected_row["leave_one_plate_out_auc"]:.3f}**

完整数字见 `metrics/strategy_summary.csv`，最终矩阵见 `matrices/final_normalized_study_matrix.tsv.gz`。
"""
    (OUT / "README.md").write_text(md, encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
