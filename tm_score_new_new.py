def process_complex(complex_dir: str, ben_scores: dict, variant_name: str):
    import os
    from glob import glob

    pdb_id = os.path.basename(complex_dir)
    print(f"\n🔍 Processing {pdb_id} [{variant_name}]")

    # Paths
    if variant_name == "combfold_high":
        combfold_dir = os.path.join(complex_dir, "combfold", "results_high", "assembled_results")
        score_file = os.path.join(complex_dir, "combfold", "results_high", "tm_score.txt")
    else:
        variant_dir = os.path.join(complex_dir, variant_name)
        combfold_dir = os.path.join(variant_dir, "results", "assembled_results")
        score_file = os.path.join(variant_dir, "tm_score.txt")

    # Ref structure
    ref_pdb = download_pdb_complex(pdb_id, complex_dir)

    # === CASE 1: If TM-score already computed, just read and return ===
    if os.path.exists(score_file):
        print(f"📄 Using existing TM-score file: {score_file}")
        best_model, best_tm, best_rmsd, ben_best_model, ben_best_tm = None, -1, -1, None, -1
        with open(score_file) as f:
            for line in f:
                if line.lower().startswith("our best model:"):
                    best_model = line.split(":")[1].strip()
                elif line.lower().startswith("tm-score:"):
                    best_tm = float(line.split(":")[1].strip())
                elif line.lower().startswith("rmsd:"):
                    best_rmsd = float(line.split(":")[1].strip())
                elif line.lower().startswith("ben's best model:"):
                    ben_best_model = line.split(":")[1].strip()
                elif "ben" in line.lower() and "tm-score" in line.lower():
                    ben_best_tm = float(line.split(":")[1].strip())
        return {
            "pdb_id": pdb_id.lower(),
            "scores": {
                variant_name: {"tm_score": best_tm, "rmsd": best_rmsd},
                "ben": {"tm_score": ben_best_tm} if ben_best_tm != -1 else {}
            }
        }

    # === CASE 2: Compute TM-score now ===
    if variant_name == "combfold_trivial":
        pattern = "output_clustered_0.pdb"
    else:
        pattern = "cb_*_output_0.pdb"

    clustered_models = sorted(glob(os.path.join(combfold_dir, pattern)))

    # Fallback for non-trivial variants if cb_* not found
    if not clustered_models and variant_name != "combfold_trivial":
        fallback_pattern = "output_clustered_0.pdb"
        clustered_models = sorted(glob(os.path.join(combfold_dir, fallback_pattern)))
        if clustered_models:
            print("✅ CombFold succeeded with output_clustered_0.pdb")

    if not clustered_models:
        print(f"⚠️  No clustered models found for {pdb_id} [{variant_name}]")
        return None

    # Align and compute TM-score
    best_tm, best_rmsd, best_model = -1, float("inf"), None
    for model in clustered_models:
        try:
            tm, rmsd = get_tm_score_rmsd_mmalign(ref_pdb, model)
            print(f"  → {os.path.basename(model)}: TM-score={tm:.5f}, RMSD={rmsd:.2f}")
            if tm > best_tm:
                best_tm = tm
                best_rmsd = rmsd
                best_model = os.path.basename(model)
        except Exception as e:
            print(f"  ❌ Error processing {model}: {e}")

    # Ben's scores
    pdb_id_upper = pdb_id.upper()
    ben_score_entries = ben_scores.get(pdb_id_upper, {}).get("scores", {})
    ben_best_tm = -1
    ben_best_model = None
    for model_name, score in ben_score_entries.items():
        tm_score = score.get("tm_score", -1)
        if tm_score > ben_best_tm:
            ben_best_tm = tm_score
            ben_best_model = model_name

    print(f"✅ Best of ours: {best_model} (TM={best_tm:.5f}, RMSD={best_rmsd:.2f})")
    if ben_best_model:
        print(f"✅ Best of Ben: {ben_best_model} (TM={ben_best_tm:.5f})")
    else:
        print("⚠️  No Ben results found for this PDB")

    # Write TM-score file
    with open(score_file, "w") as f:
        f.write(f"Our best model: {best_model}\n")
        f.write(f"TM-score: {best_tm:.5f}\n")
        f.write(f"RMSD: {best_rmsd:.2f}\n\n")
        if ben_best_model:
            f.write(f"Ben's best model: {ben_best_model}\n")
            f.write(f"TM-score: {ben_best_tm:.5f}\n")
        else:
            f.write("Ben's best model: Not found\n")

    # Return scores
    result = {
        "pdb_id": pdb_id.lower(),
        "scores": {
            variant_name: {"tm_score": best_tm, "rmsd": best_rmsd},
        }
    }
    if ben_best_tm != -1:
        result["scores"]["ben"] = {"tm_score": ben_best_tm}
    return result
