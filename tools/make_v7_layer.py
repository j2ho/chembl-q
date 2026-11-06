#!/usr/bin/env python
"""Build curated_v7 as a writable layer over curated_v6.

Only the decoy stages re-run for v7, so the 49 GB of structures is symlinked
and everything a stage writes is a real file. The assertion is the point of
this script: it refuses BEFORE any stage runs if a writable path would resolve
into the v6 release. Discovering that while stage 6 is running means v6 has
already been overwritten.

    python tools/make_v7_layer.py --src curated_v6 --dst curated_v7 \
        --weak-from /data/galaxy4/user/j2ho/tmp/v7_stage1/out
    python tools/make_v7_layer.py ... --check-only
"""

import argparse
import filecmp
import os
import shutil
import sys
from pathlib import Path

# Directories no decoy stage writes into. Symlinked, so the layer stays small.
LINK_DIRS = ("aligned", "comps", "pdb")

# Everything stage 4, 6 or 8 may rewrite. Each must be a real file or absent.
WRITABLE_PER_TARGET = (
    "decoys.tsv",            # stage 6, rewritten in place
    "actives_clustered.tsv",  # read by stage 4; copied so v7 is self-contained
    "actives.tsv",
    "weak.tsv",
    "inactives.tsv",
    "measured.tsv",
    "conflicts.tsv",
    "pocket_info.csv",
    "sequence.fasta",
)
WRITABLE_TOP = (
    "compound_pool.pkl",      # stage 4
    "chembl_targets.tsv",     # stage 8
    "best_structure.tsv",
    "passed_targets.txt",
    "pairwise_seqid.tsv",
    "pairwise_pocket_hungarian.tsv",
)


def targets_of(root: Path):
    return sorted(p.parent.name for p in root.glob("*/actives_clustered.tsv"))


def build(src: Path, dst: Path, weak_from: Path | None) -> int:
    names = targets_of(src)
    if not names:
        print(f"no targets under {src}", file=sys.stderr)
        return 1
    dst.mkdir(parents=True, exist_ok=True)

    for name in WRITABLE_TOP:
        source = src / name
        if source.exists() and not (dst / name).exists():
            shutil.copy2(source, dst / name)
    for extra in src.glob("*.tsv"):
        if not (dst / extra.name).exists():
            shutil.copy2(extra, dst / extra.name)

    n_weak = n_missing_weak = 0
    for name in names:
        s_dir, d_dir = src / name, dst / name
        d_dir.mkdir(exist_ok=True)
        for item in s_dir.iterdir():
            target = d_dir / item.name
            if item.name in LINK_DIRS:
                if not target.exists(follow_symlinks=False):
                    target.symlink_to(item.resolve())
            elif item.is_file():
                if not target.exists():
                    shutil.copy2(item, target)
        if weak_from is not None:
            weak = weak_from / name / "weak.tsv"
            if weak.exists():
                shutil.copy2(weak, d_dir / "weak.tsv")
                n_weak += 1
            else:
                n_missing_weak += 1
                # An absent weak.tsv silently narrows the reference set, so
                # write the header rather than leaving the stage to guess.
                (d_dir / "weak.tsv").write_text("chembl_id\tevidence\tsmiles\n")

    print(f"layer built: {len(names)} targets, weak.tsv from stage 1 for "
          f"{n_weak}, empty for {n_missing_weak}")
    return 0


def check(src: Path, dst: Path) -> int:
    """Every writable path is a real file; every link points outside the layer."""
    problems = []
    for name in WRITABLE_TOP:
        path = dst / name
        if path.is_symlink():
            problems.append(f"{path} is a symlink -> {os.readlink(path)}")

    names = targets_of(dst)
    if not names:
        problems.append(f"no targets under {dst}")

    for name in names:
        d_dir = dst / name
        for fname in WRITABLE_PER_TARGET:
            path = d_dir / fname
            if path.is_symlink():
                problems.append(f"{path} is a symlink -> {os.readlink(path)}")
        for link in LINK_DIRS:
            path = d_dir / link
            if path.exists(follow_symlinks=False) and not path.is_symlink():
                problems.append(f"{path} is a real directory, expected a link")

    # The decisive test: writing through the layer must not reach the source.
    probe = dst / names[0] / "decoys.tsv"
    mirror = src / names[0] / "decoys.tsv"
    if probe.exists() and mirror.exists():
        if probe.resolve() == mirror.resolve():
            problems.append(
                f"{probe} and {mirror} resolve to the same inode; writing v7 "
                f"would overwrite the v6 release")
        elif filecmp.cmp(probe, mirror, shallow=False):
            print(f"  {probe.name} is a distinct copy with identical contents, "
                  f"as expected before stage 6")

    if problems:
        print(f"REFUSING: {len(problems)} problem(s)", file=sys.stderr)
        for line in problems[:40]:
            print(f"  {line}", file=sys.stderr)
        return 1

    n_links = sum(1 for name in names for link in LINK_DIRS
                  if (dst / name / link).is_symlink())
    print(f"check passed: {len(names)} targets, {n_links} symlinked dirs, "
          f"no writable path resolves into {src}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", type=Path, required=True)
    parser.add_argument("--dst", type=Path, required=True)
    parser.add_argument("--weak-from", type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()

    if not args.check_only:
        rc = build(args.src, args.dst, args.weak_from)
        if rc:
            return rc
    return check(args.src, args.dst)


if __name__ == "__main__":
    sys.exit(main())
