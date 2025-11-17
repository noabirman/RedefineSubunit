import requests
import time

def get_phospho_sites(uniprot_id):
    """
    Fetch phosphorylation sites from UniProt JSON API
    Returns a list of dicts: {"position": int, "type": str}
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    r = requests.get(url)
    if r.status_code != 200:
        print(f"Error fetching {uniprot_id}: HTTP {r.status_code}")
        return []

    data = r.json()
    phospho_sites = []

    for feature in data.get("features", []):
        if feature.get("category") == "PTM" and "Phospho" in feature.get("description", ""):
            # Extract position safely
            loc = feature.get("location", {})
            pos = None
            if "start" in loc and "value" in loc["start"]:
                pos = loc["start"]["value"]
            elif "position" in loc and "value" in loc["position"]:
                pos = loc["position"]["value"]
            if pos is not None:
                phospho_sites.append({"position": pos, "type": feature["description"]})

    return phospho_sites

# Test
for uid in ["O75791", "P41240"]:  # replace with your UniProt IDs
    sites = get_phospho_sites(uid)
    print(uid, "->", sites)
    time.sleep(0.3)  # avoid rate limiting
