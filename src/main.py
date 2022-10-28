from warnings import catch_warnings, simplefilter
import MDAnalysis as mda
import argparse

from . import cmd_line
from . import inout
from . import apples2apples
from . import modelling


def main_align(args: argparse.Namespace):

    input_data = args.input_file
    temp = args.temp_directory
    not_aligned_sel = args.not_aligned_sel

    seqs = apples2apples.aligned_sequences(
        input_data["id"], input_data["record"], temp)

    indeces = [i-1 for i in input_data["first_residue_index"]]

    sels = apples2apples.find_apples2apples(seqs, indeces, not_aligned_sel)

    common_sels = []
    for sel, subset in zip(sels, input_data["subset"]):
        common_sel = subset.select_atoms(sel)
        common_sels.append(common_sel)

    for common_sel, ndx_file in zip(common_sels, input_data["output_ndx"]):
        with mda.selections.gromacs.SelectionWriter(ndx_file, mode='w') as ndx:
            ndx.write(common_sel, name='apples2apples')

    if input_data["output_pdb"] is not None:
        for common_sel, pdb in zip(common_sels, input_data["output_pdb"]):
            common_sel.write(pdb)


def main_model(args: argparse.Namespace):

    input_data = args.input_file
    universes = input_data["univ"]
    trajs = input_data["input_xtc"]
    print("Loading trajectories into universes")
    for pdb in universes:
        universes[pdb].load_new(trajs[pdb])

    ndxs = {pdb: [] for pdb in universes}
    print("Making selections")
    for pdb, ndx in zip(input_data["input_pdb"], input_data["output_ndx"]):
        ndxs[pdb] += inout.read_ndx(ndx)["apples2apples"]

    print("Loading trajectories to memory and translating and rotating for minimum RMSD fit")
    sels, coords, rmsds = modelling.load_to_memory(
        universes, ndxs, args.ref_frame, ref_key=input_data["input_pdb"][0])

    print()
    for key in rmsds:
        print(("For system in %s:\n"
              "rmsd min=%.4f, max=%.4f\n"
               "mean=%.4f, std=%.4f\n") % (
            key, rmsds[key].min(), rmsds[key].max(),
            rmsds[key].mean(), rmsds[key].std()))

    if (input_data["output_traj"]):
        print("Writing trajectories to file")
        modelling.write_to_files(sels, coords, input_data["output_traj"])


def main():
    parser, subparsers = cmd_line.create_main_parser()

    cmd_line.create_input_syntax_subparser(subparsers)
    cmd_line.create_align_subparser(subparsers, func=main_align)
    cmd_line.create_model_subparser(subparsers, func=main_model)

    args = parser.parse_args()

    if (not args.verbose):
        with catch_warnings():
            simplefilter("ignore", category=(DeprecationWarning, UserWarning))
            args.func(args)
    else:
        args.func(args)
