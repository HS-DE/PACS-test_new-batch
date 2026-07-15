# PACS multi-strategy benchmark

Generated: 2026-07-15 01:32:42.598627

## Data used

- Study samples after README outlier removal: 185
- Internal-QC samples: 9
- Plate-QC samples: 9
- Protein groups: 5870
- Complete HELP set: 47
- Missing-tolerant eligible HELP set: 50

## Strategies

- S00_raw: Raw log2 protein data
- S01_legacy_group2: Legacy-like independent Run/Plate plus group2 HELP and group2 median
- S02_run_only: Internal HELP Run correction only
- S03_minimal_run_plate: Run plus QC endogenous-protein Plate correction
- S04_sequential_complete: Run + Plate + complete-HELP sample correction
- S05_sequential_complete_group_median: S04 plus group2 median normalization
- S06_group_help: Run + Plate + group2-specific HELP correction
- S07_group_help_group_median: S06 plus group2 median normalization
- S08_available_help: Missing-tolerant weighted HELP correction
- S09_dual_anchor: Dual Internal HELP/protein Run anchor + Plate + available HELP
- S10_batch_model: Protein-specific protected group2 + Run + Plate model
- S11_batch_help_pc_model: Protein-specific protected model with HELP PCs
- S12_conditional_model: Apply model correction only when technical partial R2 >= 0.05

## Screening result

The screening score is not a biological truth criterion. It rewards lower residual Run/Plate association, better QC behaviour, smaller unnecessary shifts, and preservation of the raw group effect.

1. **S02_run_only** — score 32.5; Run R2=0.01845; Plate R2=0.0372; Group R2=0.01129
2. **S03_minimal_run_plate** — score 32.5; Run R2=0.01845; Plate R2=0.0372; Group R2=0.01129
3. **S06_group_help** — score 32.5; Run R2=0.01845; Plate R2=0.0372; Group R2=0.01129
4. **S00_raw** — score 35; Run R2=0.01791; Plate R2=0.0372; Group R2=0.01129
5. **S05_sequential_complete_group_median** — score 41; Run R2=0.01134; Plate R2=0.0351; Group R2=0.02172

## Interpretation rules

- Do not select a method from PCA appearance alone.
- Prefer methods that reduce Run/Plate effects in both study proteins and technical controls.
- Treat group2-dependent preprocessing strategies as sensitivity analyses unless independently validated.
- Review extreme sample factors and HELP detection counts before accepting sample-level correction.
- The final choice should combine technical-control metrics, factor stability and biological plausibility.
