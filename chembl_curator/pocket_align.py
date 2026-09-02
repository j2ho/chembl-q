# chembl_curator/pocket_align.py

"""Order-free pocket superposition: Hungarian matching plus iterative Kabsch.

Stage 5's original pocket RMSD rides on a *global* TM-align superposition and
keeps only the residue pairs TM-align matched, which is a sequence-ordered
alignment. Two pockets built from the same residue types in a different order,
or sitting in unrelated folds, cannot be matched that way even when their
geometry is nearly identical. Cross-fold pocket similarity is exactly the
leakage decoy selection cares about, so the matching has to be order-free.

Residues are paired by the Hungarian algorithm on a cost that combines a local
distance fingerprint with an amino-acid class mismatch penalty, then the pairing
and the superposition are refined against each other until the assignment stops
changing. The worst-matching tail is trimmed each round so one badly placed
residue cannot drag the fit.

Pure numpy/scipy, no structure I/O, so the geometry can be tested on synthetic
point sets.

Adapted from the cross-receptor clustering script by Su Kim
(MotifScreen_data/script/crossrec_folddsico.py). The FoldDisco index step of
that script is deliberately not reproduced: it is a candidate filter worth
having across tens of thousands of pockets, while an all-vs-all sweep of this
project's ~1.3k targets costs well under a core-hour.
"""

import math
from typing import List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

# Coarse amino-acid classes: p(ositive), n(egative), h(ydrophilic),
# b(ulk/apolar), a(romatic-tryptophan). Pairing residues across classes is
# penalised rather than forbidden.
AA_CLASS = {
    "ARG": "p", "HIS": "p", "LYS": "p", "ASP": "n", "GLU": "n",
    "ASN": "h", "GLN": "h", "SER": "h", "THR": "h", "TYR": "h",
    "ALA": "b", "CYS": "b", "GLY": "b", "ILE": "b", "LEU": "b",
    "MET": "b", "PHE": "b", "PRO": "b", "VAL": "b", "TRP": "a",
    "MSE": "b",
}

MIN_MATCHED_RESIDUES = 5
TRIM_FRACTION = 0.15
MISMATCH_PENALTY = 3.0
MAX_ITERATIONS = 20
FINGERPRINT_NEIGHBOURS = 8


def fingerprint(xyz: np.ndarray, neighbours: int = FINGERPRINT_NEIGHBOURS) -> np.ndarray:
    """Sorted distances to the nearest `neighbours` residues, per residue.

    A rotation-invariant local descriptor, so the first Hungarian round has a
    sensible cost matrix before any superposition exists.
    """
    dist = np.linalg.norm(xyz[:, None] - xyz[None, :], axis=2)
    dist.sort(axis=1)
    value = dist[:, 1:min(len(xyz), neighbours + 1)]
    if value.shape[1] < neighbours:
        pad = neighbours - value.shape[1]
        value = np.c_[value, np.repeat(value[:, -1:], pad, axis=1)]
    return value


def kabsch(mobile: np.ndarray, reference: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Rotation and translation putting `mobile` onto `reference`."""
    cm, cr = mobile.mean(0), reference.mean(0)
    u, _, vt = np.linalg.svd((mobile - cm).T @ (reference - cr))
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:          # reject reflections
        vt[-1] *= -1
        rotation = vt.T @ u.T
    return rotation, cr - cm @ rotation.T


def _cost(fp_a: np.ndarray, fp_b: np.ndarray, penalty: np.ndarray) -> np.ndarray:
    return np.linalg.norm(fp_a[:, None] - fp_b[None, :], axis=2) + penalty


def aligned_rmsd(
    pocket_a: Sequence[Tuple[str, np.ndarray]],
    pocket_b: Sequence[Tuple[str, np.ndarray]],
    min_matched: int = MIN_MATCHED_RESIDUES,
    trim_fraction: float = TRIM_FRACTION,
    mismatch_penalty: float = MISMATCH_PENALTY,
    max_iterations: int = MAX_ITERATIONS,
) -> Tuple[int, float]:
    """Superpose two pockets without using residue order.

    Each pocket is a sequence of (residue_name, CA coordinate).
    Returns (n_matched, rmsd). Raises ValueError if either pocket is empty.
    """
    if len(pocket_a) == 0 or len(pocket_b) == 0:
        raise ValueError("empty pocket")

    cls_a = np.array([AA_CLASS.get(name, "X") for name, _ in pocket_a])
    cls_b = np.array([AA_CLASS.get(name, "X") for name, _ in pocket_b])
    xyz_a = np.stack([c for _, c in pocket_a])
    xyz_b = np.stack([c for _, c in pocket_b])

    penalty = (cls_a[:, None] != cls_b) * mismatch_penalty
    floor = min(min_matched, len(xyz_a), len(xyz_b))
    cost = _cost(fingerprint(xyz_a), fingerprint(xyz_b), penalty)

    previous: Optional[Tuple] = None
    left = right = None
    rotation = np.eye(3)
    translation = np.zeros(3)

    for _ in range(max_iterations):
        left, right = linear_sum_assignment(cost)
        keep = max(floor, math.ceil(len(left) * (1.0 - trim_fraction)))
        chosen = np.argsort(cost[left, right])[:keep]
        left, right = left[chosen], right[chosen]

        rotation, translation = kabsch(xyz_b[right], xyz_a[left])
        fitted = xyz_b @ rotation.T + translation
        # Re-score against the current superposition, not the fingerprints.
        cost = np.linalg.norm(xyz_a[:, None] - fitted[None, :], axis=2) + penalty

        assignment = tuple(sorted(zip(left.tolist(), right.tolist())))
        if assignment == previous:
            break
        previous = assignment

    fitted = xyz_b[right] @ rotation.T + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - xyz_a[left]) ** 2, axis=1))))
    return len(left), rmsd


def pocket_from_structure(
    residues: Sequence[Tuple[str, np.ndarray, np.ndarray]],
    ligand_xyz: np.ndarray,
    radius: float = 8.0,
    min_residues: int = 6,
) -> Optional[List[Tuple[str, np.ndarray]]]:
    """Select pocket residues within `radius` of any ligand heavy atom.

    residues: (residue_name, CA coordinate, all heavy-atom coordinates).
    Distance is measured from every heavy atom of the residue, not from the CA
    and not from the ligand centroid, so an elongated ligand does not lose the
    residues packed against its far end.

    Returns None when the pocket is too small to align meaningfully.
    """
    if len(ligand_xyz) == 0:
        return None

    tree = cKDTree(ligand_xyz)
    pocket = [
        (name, ca)
        for name, ca, atoms in residues
        if len(atoms) and tree.query(atoms)[0].min() <= radius
    ]
    return pocket if len(pocket) >= min_residues else None
