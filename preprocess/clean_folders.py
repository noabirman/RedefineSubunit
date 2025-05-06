import os
import re
import sys
import shutil
from pathlib import Path
from datetime import datetime

def main(input_dir):
    input_dir = Path(input_dir).resolve()
    if not input_dir.is_dir():
        print(f"❌ Error: {input_dir} is not a directory.")
        return

    print(f"📂 Input directory: {input_dir}")

    # Output folder: one level above input, named "doubles"
    output_dir = input_dir.parent / "doubles"
    output_dir.mkdir(exist_ok=True)
    print(f"📁 Archive directory: {output_dir}")

    # Match folders like a_b or a_b_YYYYMMDD_HHMMSS
    folder_pattern = re.compile(r"^([a-zA-Z]+_[a-zA-Z]+)(?:_(\d{8}_\d{6}))?$")

    groups = {}
    print("\n🔍 Scanning folders...")
    for entry in input_dir.iterdir():
        if entry.is_dir():
            match = folder_pattern.match(entry.name)
            if match:
                base = match.group(1)
                timestamp_str = match.group(2)
                timestamp = None
                if timestamp_str:
                    try:
                        timestamp = datetime.strptime(timestamp_str, "%Y%m%d_%H%M%S")
                    except ValueError:
                        print(f"  ⚠️ Skipping invalid timestamp: {entry.name}")
                        continue
                else:
                    # No timestamp: treat as oldest possible
                    timestamp = datetime.min

                print(f"  ✔ Matched: {entry.name} (base='{base}', timestamp='{timestamp}')")
                groups.setdefault(base, []).append((entry, timestamp))
            else:
                print(f"  ⚠️ Skipping non-matching folder: {entry.name}")

    print("\n📊 Processing groups to keep newest folder only...")
    for base, items in groups.items():
        print(f"\n➡️ Group: {base}")
        for entry, ts in items:
            print(f"   - {entry.name} (timestamp={ts})")

        # Sort by timestamp, newest last
        items.sort(key=lambda x: x[1])
        to_keep = items[-1][0]
        to_move = [item[0] for item in items if item[0] != to_keep]

        print(f"   ✅ Keeping: {to_keep.name}")
        for folder in to_move:
            print(f"   🗂️ Moving: {folder.name} → {output_dir}")
            shutil.move(str(folder), output_dir)

    print("\n✅ Done.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python clean_folders.py /path/to/input")
        sys.exit(1)
    main(sys.argv[1])
