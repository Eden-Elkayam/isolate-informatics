"""
merge_gbk.py

For each subdirectory in a given directory, finds all .gbk/.gbff files
and merges them into a single file: <folder_name>.merged.gbk

Usage:
    python merge_gbk.py <parent_dir> [--out <output_dir>] [--ext gbk]
"""

import argparse
from pathlib import Path


def merge_gbk_files(folder, out_dir, ext):
    files = sorted(folder.glob(f"*.{ext}"))
    if not files:
        return None

    out_path = out_dir / f"{folder.name}.merged.{ext}"
    with open(out_path, 'w') as out:
        for f in files:
            text = f.read_text()
            out.write(text)
            if not text.endswith('\n'):
                out.write('\n')

    return out_path, len(files)


def main():
    parser = argparse.ArgumentParser(description="Merge GBK files per subfolder.")
    parser.add_argument("parent_dir", help="Directory containing per-genome subfolders")
    parser.add_argument("--out", help="Output directory (default: same as parent_dir)", default=None)
    parser.add_argument("--ext", help="File extension to look for (default: gbk)", default="gbk")
    args = parser.parse_args()

    parent = Path(args.parent_dir)
    out_dir = Path(args.out) if args.out else parent
    out_dir.mkdir(parents=True, exist_ok=True)

    folders = sorted([f for f in parent.iterdir() if f.is_dir()])
    if not folders:
        print(f"No subdirectories found in {parent}")
        return

    print(f"Found {len(folders)} folders in {parent}\n")
    merged_count = 0

    for folder in folders:
        result = merge_gbk_files(folder, out_dir, args.ext)
        if result:
            out_path, n = result
            print(f"  {folder.name}: merged {n} file(s) -> {out_path.name}")
            merged_count += 1
        else:
            print(f"  {folder.name}: no .{args.ext} files found, skipping")

    print(f"\nDone. {merged_count}/{len(folders)} folders merged.")


if __name__ == "__main__":
    main()
