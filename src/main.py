from warnings import catch_warnings, simplefilter
import MDAnalysis as mda

from . import inout
from . import apples2apples


def main():
    parser = inout.create_parser()

    args = parser.parse_args()

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
