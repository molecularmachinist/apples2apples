from typing import Dict, List, Tuple, Union

import pathlib
import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import align

from . import utils


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


def fit_and_write(univs: Dict[str, mda.Universe],
                  ndxs: Dict[str, List[int]],
                  outtrajs: Dict[str, str],
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
    for j, key in enumerate(univs):
        u = univs[key]
        sel = sels[key]
        rmsd[key] = np.empty(len(u.trajectory))
        len_trj = len(u.trajectory)
        num_sys = len(univs)
        with mda.Writer(outtrajs[key], n_atoms=sel.n_atoms) as w:
            for i, ts in enumerate(u.trajectory):
                sel.positions, rmsd[key][i] = fit_pos(sel.positions, ref_pos)
                if (i == 0):
                    sel.write(str(
                        pathlib.Path(outtrajs[key]).with_suffix(".pdb")))
                    sel.write(str(
                        pathlib.Path(outtrajs[key]).with_suffix(".ndx")),
                        name='apples2apples')
                print(f"trajectory {j+1}/{num_sys}: {i+1}/{len_trj}", end="\r")
                w.write(sel)
            print("")

    return rmsd


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
