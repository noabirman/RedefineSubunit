#!/usr/bin/env python3
import json
import requests
import argparse
import os
import time

# -----------------------------
# Function to extract canonical UniProt accession
# -----------------------------
def get_accession(orig_id):
    # Handles formats like sp|P12345|PROT_HUMAN
    if "|" in orig_id:
        return orig_id.split("|")[1]
    return orig_id

# -----------------------------
# Function to fetch phospho sites from UniProt JSON API
# -----------------------------
def get_phospho_sites(uniprot_id):
    accession = get_accession(uniprot_id)
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    r = requests.get(url)
    if r.status_code != 200:
        print(f"Warning: UniProt request failed for {uniprot_id} (accession {accession}) - status {r.status_code}")
        return []

    data = r.json()
    phospho_sites = []

    for feature in data.get("features", []):
        # Only PTM features
        if feature.get("category") == "PTM" and "Phospho" in feature.get("description", ""):
            # Extract position safely
            pos = None
            loc = feature.get("location", {})
            if "start" in loc and "value" in loc["start"]:
                pos = loc["start"]["value"]
            elif "position" in loc and "value" in loc["position"]:
                pos = loc["position"]["value"]
            if pos is not None:
                phospho_sites.append({"position": pos, "type": feature["description"]})

    return phospho_sites

# -----------------------------
# Main script
# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map phosphorylation sites to subunits")
    parser.add_argument('--subunits_info', required=True, help='Path to subunits_info.json')
    parser.add_argument('--chain_map', required=True, help='Path to chain_id_mapping.json')
    parser.add_argument('--output', default='subunit_phospho.json', help='Output JSON file path')
    args = parser.parse_args()

    # Load input files
    if not os.path.exists(args.subunits_info):
        raise FileNotFoundError(f"{args.subunits_info} not found")
    if not os.path.exists(args.chain_map):
        raise FileNotFoundError(f"{args.chain_map} not found")

    with open(args.subunits_info) as f:
        subunits = json.load(f)

    with open(args.chain_map) as f:
        chain_map = json.load(f)

    subunit_phospho = {}

    # Loop over each protein
    for orig_id, info in subunits.items():
        accession = get_accession(orig_id)
        print(f"Processing {orig_id} (accession {accession})...")
        phospho_sites = get_phospho_sites(orig_id)
        print(f"  Found {len(phospho_sites)} phospho sites in UniProt")
        time.sleep(0.3)  # rate limiting

        for site in phospho_sites:
            site_pos = site["position"]

            for subunit_key, subunit_info in chain_map.items():
                # Only consider subunits that originate from this protein
                if subunit_info["chain_id"] not in info["chain_names"]:
                    continue

                start, end = subunit_info["start"], subunit_info["end"]
                # Adjust start if 0-based
                start_corrected = start + 1

                if start_corrected <= site_pos <= end:
                    relative_pos = site_pos - start_corrected + 1
                    subunit_phospho.setdefault(subunit_key, []).append({
                        "position": relative_pos,
                        "type": site["type"]
                    })
                    print(f"  Mapping UniProt site {site_pos} -> {subunit_key} position {relative_pos}")

    # Save output
    with open(args.output, "w") as f:
        json.dump(subunit_phospho, f, indent=2)

    print(f"\nPhosphorylation sites mapped to subunits saved as {args.output}")
