import pathlib
from typing import Dict, List, Tuple, Union
import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np


def fit_pos(pos: np.ndarray, ref: np.ndarray) -> Tuple[
        np.ndarray, np.ndarray]:
    pos_c = pos.mean(axis=0)
    ref_c = ref.mean(axis=0)

    newpos = pos-pos_c
    cenref = ref-ref_c
    R, rmsd = align.rotation_matrix(newpos, cenref)
    newpos = newpos @ R.T
    newpos += ref_c
    return newpos, rmsd


def load_to_memory(univs: Dict[str, mda.Universe],
                   ndxs: Dict[str, List[int]],
                   ref_frame: int,
                   ref_key: str
                   ) -> Tuple[
        Dict[str, mda.AtomGroup],
        Dict[str, np.ndarray],
        Dict[str, np.ndarray]]:
    trjs = {}
    rmsd = {}
    sels = {key: univs[key].atoms[np.array(ndxs[key])-1] for key in univs}
    univs[ref_key].trajectory[ref_frame]
    ref_pos = sels[ref_key].positions
    for key in univs:
        u = univs[key]
        sel = sels[key]
        res = [fit_pos(sel.positions, ref_pos) for ts in u.trajectory]
        trjs[key] = np.array([pos.flatten() for pos, _ in res])
        rmsd[key] = np.array([rmsd for _, rmsd in res])
    return sels, trjs, rmsd


def write_to_files(sels: Dict[str, mda.AtomGroup],
                   trajs: Dict[str, np.ndarray],
                   outtrajs: Dict[str, str]):

    for key in trajs:
        trj = trajs[key]
        sel = sels[key]
        sel.positions = trj[0].reshape((-1, 3))
        sel.write(str(pathlib.Path(outtrajs[key]).with_suffix(".pdb")))
        sel.write(str(pathlib.Path(outtrajs[key]).with_suffix(".ndx")),
                  name='apples2apples')
        with mda.Writer(outtrajs[key], n_atoms=trj.shape[1]) as w:
            for arr in trj:
                sel.positions = arr.reshape(-1, 3)
                w.write(sels[key])
