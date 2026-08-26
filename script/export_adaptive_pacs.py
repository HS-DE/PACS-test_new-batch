#!/usr/bin/env python3
import os
import json
import gzip
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "adaptive_export")
os.makedirs(OUT, exist_ok=True)

OUTLIERS = [
    "WJ-1-5","WJ-1-7","WJ-1-19","WJ-1-21","WJ-1-44","WJ-1-67",
    "WJ-1-90","WJ-1-101","WJ-1-102","WJ-1-103","WJ-1-104",
    "WJ-1-105","WJ-1-112","WJ-1-146","WJ-1-148","WJ-1-157",
    "WJ-1-171","WJ-1-176","WJ-1-177","WJ-1-178","WJ-1-187",
    "WJ-1-203","WJ-1-208","WJ-1-209","WJ-1-212","WJ-1-213",
    "WJ-1-43","WJ-1-59","WJ-1-83"
]

PRIMARY_HELP = [
    "IQNILTEEPK",
    "AGALNSNDAFVLK",
    "LAPLAEDVR",
    "GSESGIFTNTK",
    "NIQSLEVIGK",
]

# Frozen development parameters from the adaptive PACS optimization.
INSTRUMENT_STRENGTH = 0.50
PLATE_Z_THRESHOLD = 0.50
PLATE_STRENGTH = 0.25
SAMPLE_HELP_Z_THRESHOLD = 1.75
SAMPLE_HELP_MAX_ALPHA = 0.35
SAMPLE_HELP_COHERENCE_MIN = 0.70
PROTEIN_SHIFT_Z_MIN = 0.50


def robust_mad(x, scale=False):
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan
    med = np.median(a)
    mad = np.median(np.abs(a - med))
    return mad * (1.4826 if scale else 1.0)


def robust_z(s):
    s = pd.Series(s, dtype=float)
    med = np.nanmedian(s.values)
    mad = robust_mad(s.values, scale=True)
    if not np.isfinite(mad) or mad <= 1e-12:
        return pd.Series(np.zeros(len(s)), index=s.index, dtype=float)
    return (s - med) / mad


def classify_qc_role(row):
    txt = " | ".join(str(v).lower() for v in row.values)
    if "blank" in txt:
        return "blank"
    if "neat" in txt or "qc2" in txt:
        return "neat"
    if "internal" in txt or "calibrator" in txt or "qc3" in txt:
        return "internal_qc"
    if "qc1" in txt or "plate qc" in txt or "plate_qc" in txt:
        return "plate_qc"
    return "other_control"


def prepare_metadata():
    sm = pd.read_excel(os.path.join(DATA, "metadata", "all_sample_metadata.xlsx"))
    qm = pd.read_excel(os.path.join(DATA, "metadata", "all_QC_metadata.xlsx"))
    # Preserve the repository's canonical columns.
    sm["sample_id"] = sm["sample_id"].astype(str).str.strip()
    qm["sample_id"] = qm["sample_id"].astype(str).str.strip()
    sm = sm[~sm["sample_id"].isin(OUTLIERS)].copy()
    qm["sample_role"] = qm.apply(classify_qc_role, axis=1)
    return sm, qm


def load_protein():
    pg = pd.read_csv(os.path.join(DATA, "pg_matrix", "pg_matrix.tsv"), sep="\t")
    anno = pg.iloc[:, :5].copy()
    mat = pg.iloc[:, 5:].apply(pd.to_numeric, errors="coerce")
    mat.index = pg.iloc[:, 0].astype(str)
    # Same treatment as repo: non-positive/non-finite -> missing, then log2.
    mat = mat.where(np.isfinite(mat) & (mat > 0))
    log2mat = np.log2(mat)
    return pg, anno, mat, log2mat


def load_help(known_ids):
    h = pd.read_excel(os.path.join(DATA, "HELP", "ALL-HELP.xlsx"))
    # In this repo HELP rows are peptides and columns are samples.
    if sum(c in known_ids for c in h.columns) >= sum(str(v) in known_ids for v in h.iloc[:, 0]):
        pep = h.iloc[:, 0].astype(str)
        hm = h.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        hm.index = pep
    else:
        ids = h.iloc[:, 0].astype(str)
        hm = h.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").T
        hm.columns = ids
    hm = hm.where(np.isfinite(hm) & (hm > 0))
    med = np.nanmedian(hm.values)
    if med > 100:
        hm = np.log2(hm)
    return hm


def panel_factors(H, panel, ids, ref_ids):
    panel = [p for p in panel if p in H.index]
    ids = [s for s in ids if s in H.columns]
    ref_ids = [s for s in ref_ids if s in H.columns]
    sub = H.loc[panel, ids]
    ref = H.loc[panel, ref_ids].median(axis=1, skipna=True)
    dev = sub.sub(ref, axis=0)
    fac = dev.median(axis=0, skipna=True)
    coherence = pd.Series(index=ids, dtype=float)
    n_help = pd.Series(index=ids, dtype=float)
    for s in ids:
        d = dev[s].dropna()
        n_help.loc[s] = len(d)
        if len(d) == 0 or not np.isfinite(fac.loc[s]) or fac.loc[s] == 0:
            coherence.loc[s] = 0.0 if len(d) == 0 else 1.0
        else:
            coherence.loc[s] = np.mean(np.sign(d.values) == np.sign(fac.loc[s]))
    return fac, coherence, n_help, dev


def protein_plate_offsets(qc_log, qc_meta):
    """Estimate protein-specific pooled-QC plate offsets and their technical SNR."""
    plates = sorted(pd.unique(qc_meta["Plate"]))
    global_ref = qc_log.median(axis=1, skipna=True)
    offsets = {}
    zscores = {}
    within = pd.DataFrame(index=qc_log.index)

    # Robust within-plate noise is estimated from all QC residuals around each plate median.
    residual_cols = []
    for plate in plates:
        ids = qc_meta.loc[qc_meta["Plate"] == plate, "sample_id"].tolist()
        ids = [s for s in ids if s in qc_log.columns]
        if not ids:
            continue
        pmed = qc_log[ids].median(axis=1, skipna=True)
        offsets[plate] = pmed - global_ref
        residual_cols.append(qc_log[ids].sub(pmed, axis=0))
    if residual_cols:
        residuals = pd.concat(residual_cols, axis=1)
        within_mad = residuals.apply(lambda r: robust_mad(r.values, scale=True), axis=1)
    else:
        within_mad = pd.Series(np.nan, index=qc_log.index)

    for plate in plates:
        if plate not in offsets:
            continue
        n = int((qc_meta["Plate"] == plate).sum())
        se = (within_mad / np.sqrt(max(n, 1))).clip(lower=0.05)
        zscores[plate] = offsets[plate].abs() / se
    return offsets, zscores, within_mad


def apply_protein_plate_corr(X, study_meta, qc_log, qc_meta, zthr=0.5, strength=0.25):
    Y = X.copy()
    offsets, zscores, within_mad = protein_plate_offsets(qc_log, qc_meta)
    applied_count = pd.Series(0, index=X.index, dtype=int)
    for plate, off in offsets.items():
        z = zscores[plate]
        use = z >= zthr
        ids = study_meta.loc[study_meta["Plate"] == plate, "sample_id"].tolist()
        ids = [s for s in ids if s in Y.columns]
        if ids:
            correction = (strength * off.where(use, 0.0)).fillna(0.0)
            Y.loc[:, ids] = Y.loc[:, ids].sub(correction, axis=0)
            applied_count += use.fillna(False).astype(int)
    return Y, offsets, zscores, within_mad, applied_count


def protein_factor_for(X):
    # Label-free sample-wide endogenous shift relative to each protein's study median.
    ref = X.median(axis=1, skipna=True)
    dev = X.sub(ref, axis=0)
    return dev.median(axis=0, skipna=True)


def row_mad_df(df):
    arr = df.to_numpy(dtype=float)
    meds = np.nanmedian(arr, axis=1)
    return np.nanmedian(np.abs(arr - meds[:, None]), axis=1)


def geometry_metrics(X, meta, feature_idx):
    from numpy.linalg import svd
    M = X.iloc[feature_idx].T.copy()
    # training-free median imputation for descriptive geometry only
    M = M.apply(lambda c: c.fillna(c.median()), axis=0)
    M = M.loc[:, M.std(axis=0, ddof=0) > 0]
    A = M.to_numpy(dtype=float)
    A = (A - np.mean(A, axis=0)) / np.std(A, axis=0, ddof=0)
    U, S, Vt = svd(A, full_matrices=False)
    k = min(10, U.shape[1])
    Z = U[:, :k] * S[:k]
    m = meta.set_index("sample_id").loc[M.index]
    groups = m["group2"].astype(str).values
    uniq = list(pd.unique(groups))
    centroids = {g: Z[groups == g].mean(axis=0) for g in uniq}
    within = {}
    for g in uniq:
        d = np.linalg.norm(Z[groups == g] - centroids[g], axis=1)
        within[g] = float(np.median(d))
    if len(uniq) >= 2:
        between = float(np.linalg.norm(centroids[uniq[0]] - centroids[uniq[1]]))
        mean_within = float(np.mean([within[uniq[0]], within[uniq[1]]]))
        ratio = between / mean_within if mean_within > 0 else np.nan
    else:
        between = ratio = np.nan

    def eta2_numeric(pc, labels):
        labels = np.asarray(labels).astype(str)
        grand = np.mean(pc)
        sst = np.sum((pc-grand)**2)
        if sst <= 0: return 0.0
        ssb = 0.0
        for g in pd.unique(labels):
            vals = pc[labels == g]
            ssb += len(vals) * (np.mean(vals)-grand)**2
        return ssb/sst

    # weighted average across first five PCs using explained variance weights
    var = S**2
    weights = var[:min(5,len(var))] / np.sum(var[:min(5,len(var))])
    plate_eta = sum(weights[j]*eta2_numeric(Z[:,j], m["Plate"].values) for j in range(len(weights)))
    run_eta = sum(weights[j]*eta2_numeric(Z[:,j], m["Run"].values) for j in range(len(weights)))
    return within, between, ratio, plate_eta, run_eta


def save_sample_x_protein(path, X, sample_ids):
    out = X.loc[:, sample_ids].T.copy()
    out.insert(0, "sample_id", out.index)
    out.to_csv(path, sep="\t", index=False, na_rep="NA", compression="gzip")


# ------------------------- load -------------------------
sm, qm = prepare_metadata()
pg, anno, raw_intensity, Xlog = load_protein()
known_ids = set(sm["sample_id"]).union(set(qm["sample_id"]))
H = load_help(known_ids)

study_ids = [s for s in sm["sample_id"] if s in Xlog.columns and s in H.columns]
sm = sm.set_index("sample_id").loc[study_ids].reset_index()
qm_internal = qm[qm["sample_role"] == "internal_qc"].copy()
qm_plate = qm[qm["sample_role"] == "plate_qc"].copy()
internal_ids = [s for s in qm_internal["sample_id"] if s in Xlog.columns and s in H.columns]
plate_ids = [s for s in qm_plate["sample_id"] if s in Xlog.columns and s in H.columns]
qm_internal = qm_internal.set_index("sample_id").loc[internal_ids].reset_index()
qm_plate = qm_plate.set_index("sample_id").loc[plate_ids].reset_index()

Xraw = Xlog.loc[:, study_ids].copy()
Qplate = Xlog.loc[:, plate_ids].copy()

# ------------------------- HELP factors -------------------------
fi, ci, ni, _ = panel_factors(H, PRIMARY_HELP, internal_ids, internal_ids)
run_factor = {}
for run, g in qm_internal.assign(help_factor=qm_internal["sample_id"].map(fi)).groupby("Run"):
    run_factor[run] = float(np.nanmedian(g["help_factor"]))

fq, cq, nq, _ = panel_factors(H, PRIMARY_HELP, plate_ids, plate_ids)
plate_help_factor = {}
for plate, g in qm_plate.assign(help_factor=qm_plate["sample_id"].map(fq)).groupby("Plate"):
    plate_help_factor[plate] = float(np.nanmedian(g["help_factor"]))

fs, cs, ns, _ = panel_factors(H, PRIMARY_HELP, study_ids, study_ids)
help_residual = fs.copy()
sm_index = sm.set_index("sample_id")
for sid in study_ids:
    r = sm_index.loc[sid, "Run"]
    p = sm_index.loc[sid, "Plate"]
    help_residual.loc[sid] = fs.loc[sid] - run_factor.get(r, 0.0) - plate_help_factor.get(p, 0.0)
help_z = robust_z(help_residual)

# ------------------------- protein plate shrinkage -------------------------
PX, plate_offsets, plate_z, plate_noise, plate_applied = apply_protein_plate_corr(
    Xraw, sm, Qplate, qm_plate,
    zthr=PLATE_Z_THRESHOLD,
    strength=PLATE_STRENGTH
)

# ------------------------- partial instrument correction -------------------------
rv = pd.Series({sid: run_factor.get(sm_index.loc[sid, "Run"], 0.0) for sid in study_ids})
BX = PX.sub(INSTRUMENT_STRENGTH * rv, axis=1)

# ------------------------- adaptive study-sample correction -------------------------
protein_factor = protein_factor_for(BX)
protein_z = robust_z(protein_factor)
az = help_z.abs()
alpha = pd.Series(0.0, index=study_ids)
mask = az > SAMPLE_HELP_Z_THRESHOLD
alpha.loc[mask] = SAMPLE_HELP_MAX_ALPHA * (1.0 - (SAMPLE_HELP_Z_THRESHOLD / az.loc[mask]) ** 2)

gate = (
    (cs.reindex(study_ids).fillna(0.0) >= SAMPLE_HELP_COHERENCE_MIN)
    & (np.sign(help_residual.reindex(study_ids).fillna(0.0)) == np.sign(protein_factor.reindex(study_ids).fillna(0.0)))
    & (protein_z.reindex(study_ids).abs().fillna(0.0) >= PROTEIN_SHIFT_Z_MIN)
)
alpha = alpha.where(gate, 0.0).clip(lower=0.0, upper=SAMPLE_HELP_MAX_ALPHA)
Xfinal = BX.sub(alpha * help_residual.reindex(study_ids), axis=1)

# ------------------------- label-free visualization features -------------------------
study_det = raw_intensity.loc[:, study_ids].notna().mean(axis=1).values
qc_det = raw_intensity.loc[:, plate_ids].notna().mean(axis=1).values
study_mad = row_mad_df(Xraw)
qc_mad = row_mad_df(Qplate)
ratio = study_mad / np.maximum(qc_mad, 0.03)
fmask = (study_det >= 0.70) & (qc_det >= (8/9)) & np.isfinite(ratio)
elig = np.where(fmask)[0]
order = elig[np.argsort(ratio[fmask])]
top500_idx = order[-min(500, len(order)):]
top500_names = Xraw.index[top500_idx].tolist()

# ------------------------- outputs -------------------------
save_sample_x_protein(os.path.join(OUT, "01_before_log2_sample_x_protein.tsv.gz"), Xraw, study_ids)
save_sample_x_protein(os.path.join(OUT, "02_after_adaptive_pacs_sample_x_protein.tsv.gz"), Xfinal, study_ids)
save_sample_x_protein(os.path.join(OUT, "03_before_log2_top500_sample_x_protein.tsv.gz"), Xraw.loc[top500_names], study_ids)
save_sample_x_protein(os.path.join(OUT, "04_after_adaptive_pacs_top500_sample_x_protein.tsv.gz"), Xfinal.loc[top500_names], study_ids)

meta_out = sm.copy()
meta_out.to_csv(os.path.join(OUT, "05_sample_metadata.tsv"), sep="\t", index=False)

corr = sm[["sample_id", "sample_name", "Plate", "Run", "group1", "group2"]].copy()
corr["HELP_factor_raw"] = corr["sample_id"].map(fs)
corr["HELP_factor_residual"] = corr["sample_id"].map(help_residual)
corr["HELP_z"] = corr["sample_id"].map(help_z)
corr["HELP_coherence"] = corr["sample_id"].map(cs)
corr["n_HELP"] = corr["sample_id"].map(ns)
corr["protein_global_factor"] = corr["sample_id"].map(protein_factor)
corr["protein_global_z"] = corr["sample_id"].map(protein_z)
corr["alpha"] = corr["sample_id"].map(alpha)
corr["sample_HELP_correction_applied"] = corr["alpha"] > 0
corr.to_csv(os.path.join(OUT, "06_sample_correction_factors.tsv"), sep="\t", index=False)

feature_info = anno.copy()
feature_info["study_detection"] = study_det
feature_info["plate_QC_detection"] = qc_det
feature_info["study_MAD_log2"] = study_mad
feature_info["plate_QC_MAD_log2"] = qc_mad
feature_info["biology_to_technical_ratio"] = ratio
feature_info["top500_visualization"] = feature_info.iloc[:,0].astype(str).isin(top500_names)
feature_info["plate_correction_applied_n_plates"] = plate_applied.values
feature_info.to_csv(os.path.join(OUT, "07_feature_annotation_and_qc.tsv.gz"), sep="\t", index=False, compression="gzip")

# technical factor summaries
pd.DataFrame([{"Run":k,"HELP_instrument_factor":v,"applied_strength":INSTRUMENT_STRENGTH} for k,v in run_factor.items()]).to_csv(
    os.path.join(OUT, "08_instrument_factors.tsv"), sep="\t", index=False)
pd.DataFrame([{"Plate":k,"HELP_plate_factor":v} for k,v in plate_help_factor.items()]).to_csv(
    os.path.join(OUT, "09_plate_HELP_factors.tsv"), sep="\t", index=False)

# descriptive validation metrics only; group label was not used to derive factors.
metrics_rows = []
for nfeat in [300, 500, 1000]:
    use_idx = order[-min(nfeat, len(order)):]
    for label, mat in [("before", Xraw), ("after", Xfinal)]:
        within, between, sep, pe, re = geometry_metrics(mat, sm, use_idx)
        row = {
            "feature_set_n": min(nfeat, len(order)), "matrix": label,
            "between_distance": between, "between_within_ratio": sep,
            "plate_eta2": pe, "run_eta2": re,
        }
        for g,v in within.items(): row[f"within_{g}"] = v
        metrics_rows.append(row)
metrics = pd.DataFrame(metrics_rows)
metrics.to_csv(os.path.join(OUT, "10_validation_metrics.tsv"), sep="\t", index=False)

params = {
    "primary_HELP": PRIMARY_HELP,
    "instrument_strength": INSTRUMENT_STRENGTH,
    "plate_z_threshold": PLATE_Z_THRESHOLD,
    "plate_strength": PLATE_STRENGTH,
    "sample_HELP_z_threshold": SAMPLE_HELP_Z_THRESHOLD,
    "sample_HELP_max_alpha": SAMPLE_HELP_MAX_ALPHA,
    "sample_HELP_coherence_min": SAMPLE_HELP_COHERENCE_MIN,
    "protein_shift_z_min": PROTEIN_SHIFT_Z_MIN,
    "n_study_samples": len(study_ids),
    "n_sample_HELP_corrected": int((alpha > 0).sum()),
    "matrix_scale": "log2 positive protein intensity; missing/non-positive retained as NA",
    "matrix_orientation": "sample x protein",
}
with open(os.path.join(OUT, "11_parameters.json"), "w") as f:
    json.dump(params, f, indent=2)

with open(os.path.join(OUT, "README.txt"), "w") as f:
    f.write("Adaptive PACS export\n")
    f.write("====================\n\n")
    f.write("01_before_log2_sample_x_protein.tsv.gz: pre-normalization study matrix after the predefined outlier removal; log2 scale.\n")
    f.write("02_after_adaptive_pacs_sample_x_protein.tsv.gz: final normalized matrix; log2 scale.\n")
    f.write("03/04: the same before/after matrices restricted to 500 label-free biology-to-technical features for quick visualization.\n")
    f.write("05_sample_metadata.tsv: group/Plate/Run metadata for plotting.\n")
    f.write("06_sample_correction_factors.tsv: HELP residual factor, robust Z, coherence, endogenous shift and applied alpha.\n")
    f.write("07_feature_annotation_and_qc.tsv.gz: protein annotation and QC-derived feature metrics.\n")
    f.write("10_validation_metrics.tsv: descriptive pre/post geometry and residual Plate/Run association. Labels are used only here for validation, not normalization.\n")

print("Adaptive PACS export complete")
print(json.dumps(params, indent=2))
print(metrics.to_string(index=False))
