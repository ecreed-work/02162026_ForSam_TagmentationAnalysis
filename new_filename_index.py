"""
Copy and re-index TSV files from one dataset into a new directory so that
their numeric indices continue sequentially from a previous dataset.

Behavior:
- Reads files whose names follow the pattern:
    {index}_S{index}_L001_R1_001.insertions.5p.tsv
- For each source file, creates a copy (does not modify or delete originals).
- Adds a fixed offset to the numeric index (e.g., +41) so numbering continues
  from where an earlier dataset ended.
- Writes the renamed copies to a user-specified destination folder.
- Preserves file contents and metadata.

Use case:
Combine two datasets that share identical naming conventions but restart
their numbering at 1, avoiding filename collisions while keeping originals.

Execute:
python new_filename_index.py \
  --src_dir /Users/ecreed/Desktop/02052026_TagmentatoinAnalysis_251218_260117_260203_combined/260203_R1_tsv_files \
  --dst_dir /Users/ecreed/Desktop/02052026_TagmentatoinAnalysis_251218_260117_260203_combined/reindexed_260203_R1_tsv_files \
  --offset 80

"""
import os
import shutil
import argparse


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--src_dir", required=True, help="Folder containing input TSV files")
    p.add_argument("--dst_dir", required=True, help="Folder to write renamed copies")
    p.add_argument("--offset", type=int, required=True, help="Add this to each index (e.g., 41)")
    p.add_argument("--start", type=int, default=1, help="First source index (default: 1)")
    p.add_argument("--end", type=int, default=40, help="Last source index inclusive (default: 40)")
    args = p.parse_args()

    os.makedirs(args.dst_dir, exist_ok=True)

    for i in range(args.start, args.end + 1):
        old_path = os.path.join(args.src_dir, f"{i}_S{i}_L001_R1_001.insertions.5p.tsv")
        new_i = i + args.offset
        new_path = os.path.join(args.dst_dir, f"{new_i}_S{new_i}_L001_R1_001.insertions.5p.tsv")
        shutil.copy2(old_path, new_path)


if __name__ == "__main__":
    main()

