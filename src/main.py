from warnings import catch_warnings, simplefilter
import MDAnalysis as mda
import argparse

from . import cmd_line
from . import apples2apples


def main_align(args: argparse.Namespace):

    ids, universes, records, subsets, first_residue_indexes, output_ndxs, output_pdbs = args.input_file
    temp = args.temp_directory
    not_aligned_sel = args.not_aligned_sel

    seqs = apples2apples.aligned_sequences(ids, records, temp)

    indeces = [i-1 for i in first_residue_indexes]

    sels = apples2apples.find_apples2apples(seqs, indeces, not_aligned_sel)

    common_sels = []
    for sel, subset in zip(sels, subsets):
        common_sel = subset.select_atoms(sel)
        common_sels.append(common_sel)

    if output_ndxs is not None:
        for common_sel, ndx_file in zip(common_sels, output_ndxs):
            with mda.selections.gromacs.SelectionWriter(ndx_file, mode='w') as ndx:
                ndx.write(common_sel, name='apples2apples')

    if output_pdbs is not None:
        with catch_warnings():
            simplefilter("ignore", category=(DeprecationWarning, UserWarning))

            for common_sel, pdb in zip(common_sels, output_pdbs):
                common_sel.write(pdb)


def main_model(args: argparse.Namespace):
    pass


def main():
    parser, subparsers = cmd_line.create_main_parser()

    cmd_line.create_input_syntax_subparser(subparsers, func=main_align)
    cmd_line.create_align_subparser(subparsers, func=main_align)
    cmd_line.create_model_subparser(subparsers, func=main_model)

    args = parser.parse_args()
    args.func(args)
