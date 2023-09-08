from warnings import catch_warnings, simplefilter
import argparse

from . import cmd_line
from . import inout
from . import apples2apples
from . import fitting
from . import utils


def main_align(args: argparse.Namespace):
    """
    Main function for the align command

    Parameters:
    -----------
    args: Namespace
        the argparse Namespace that has all commandline options
    """

    input_data = args.input_file
    temp = args.temp_directory
    not_aligned_sel = args.not_aligned_sel

    seqs = apples2apples.aligned_sequences(input_data["record"], temp)

    if (args.resid):
        resids = [[r.resid for r in sel.residues]
                  for sel in input_data["subset"]]
    else:
        resids = [[r.resindex for r in sel.residues]
                  for sel in input_data["subset"]]

    sels = apples2apples.find_apples2apples(
        seqs, resids, not_aligned_sel, args.resid
    )

    common_sels = []
    for i, (sel, subset) in enumerate(zip(sels, input_data["subset"])):
        common_sel = subset.select_atoms(sel)
        common_sels.append(common_sel)
        print(f"{utils.ordinal_str(i+1)} sequence: "
              f"{len(subset)} atoms -> {len(common_sel)} atoms")
        if (args.print_sels):
            print(sel)

    print("Writing output index files")
    for common_sel, ndx_file in zip(common_sels, input_data["output_ndx"]):
        if (args.make_dirs):
            ndx_file.parent.mkdir(parents=True, exist_ok=True)
        elif (not ndx_file.parent.is_dir()):
            raise FileNotFoundError(f"Parent directory of {ndx_file} does not "
                                    "exists, but \"--no-mkdir\" was given.")
        common_sel.write(str(ndx_file), name='apples2apples')

    if input_data["output_pdb"] is not None:
        print("Writing output pdbs")
        for common_sel, pdb in zip(common_sels, input_data["output_pdb"]):
            if (args.make_dirs):
                pdb.parent.mkdir(parents=True, exist_ok=True)
            elif (not pdb.parent.is_dir()):
                raise FileNotFoundError(f"Parent directory of {pdb} does not "
                                        "exists, but \"--no-mkdir\" was given.")
            common_sel.write(pdb)


def main_fit(args: argparse.Namespace):
    """
    Main function for the fit command

    Parameters:
    -----------
    args: Namespace
        the argparse Namespace that has all commandline options
    """

    input_data = args.input_file
    universes = input_data["univ"]
    trajs = input_data["input_xtc"]
    print("Loading trajectories into universes")
    for pdb in universes:
        universes[pdb].load_new(str(trajs[pdb]))

    ndxs = {pdb: [] for pdb in universes}
    print("Making selections")
    for pdb, ndx in zip(input_data["input_pdb"], input_data["output_ndx"]):
        ndxs[pdb] += inout.read_ndx(ndx)["apples2apples"]

    print("Starting to read and write trajectories, translating and rotating for minimum RMSD fit")
    rmsds = fitting.fit_and_write(univs=universes,
                                  ndxs=ndxs,
                                  outtrajs=input_data["output_traj"],
                                  ref_frame=args.ref_frame,
                                  ref_key=input_data["input_pdb"][0],
                                  make_dirs=args.make_dirs)

    print()
    for key in rmsds:
        print(("For system in %s:\n"
              "rmsd min=%.4f, max=%.4f\n"
               "mean=%.4f, std=%.4f\n") % (
            key, rmsds[key].min(), rmsds[key].max(),
            rmsds[key].mean(), rmsds[key].std()))


def main():
    """
    Main function of the program
    """
    parser, subparsers = cmd_line.create_main_parser()

    cmd_line.create_align_subparser(subparsers, func=main_align)
    cmd_line.create_fit_subparser(subparsers, func=main_fit)

    args = parser.parse_args()

    if (not args.verbose):
        with catch_warnings():
            simplefilter("ignore", category=(DeprecationWarning, UserWarning))
            args.func(args)
    else:
        args.func(args)
