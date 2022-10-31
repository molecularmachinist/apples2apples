import argparse
import sys
from typing import Tuple

from . import inout


def create_main_parser() -> Tuple[argparse.ArgumentParser, argparse._SubParsersAction]:
    """
    Create the main parser of the program and a subparser for it.
    """
    description = "apples2apples: A tool for comparing apples to oranges."
    epilog = "For more detailed explanations, see the README.md"
    parser = argparse.ArgumentParser(
        prog="apples2apples", description=description, epilog=epilog)

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Unsilence warnings such as missing attributes when writing PDB files.")

    msg = "Specify what the program should do. Only one command should be given per run. 'apples2apples <cmd> -h' to get command specific option."
    subparsers = parser.add_subparsers(description=msg)
    subparsers.required = True
    subparsers.dest = 'command'

    return parser, subparsers


def create_input_syntax_subparser(subparsers: argparse._SubParsersAction, **kwargs) -> argparse.ArgumentParser:
    """
    Create the subparser for the input_syntax command. **kwargs are added to set_defaults
    """
    syntax = '''
# Input file example.
# The symbol # is used for comment lines.
# The first non-empty line after these comments includes column titles.
# The columns are separated from eachother by pipe symbol |.
# The columns input_pdbs, first_residue_index and output_ndx
# are required. In the input_pdbs at least 2 pdb files are required.
# If selections are not given, they are assumed to be "all".
# See  https://docs.mdanalysis.org/stable/documentation_pages/selections.html
# for syntax of the selections column.
# The columns output_ndx and output_pdb are to specify output files in alignment,
# while the former is used as input in fitting.
# Note that the pipe symbols don't have to be aligned.
# The below line marks which columns are used by which command, a=align and f=fit.
# Added asterisk means its required for that command

# a* f*   | f*        |  a                       | a* f*      | a          | f*

input_pdb | input_xtc |  selection               | output_ndx | output_pdb | output_traj
1.pdb     | 1.pdb     |  segid A and resid 40:60 | 1.ndx      | 1out.pdb   | 1aligned.xtc
2.pdb     | 2.pdb     |  segid B and resid 10:55 | 2.ndx      | 2out.pdb   | 2aligned.xtc
3.pdb     | 3.pdb     |  protein                 | 3.ndx      | 3out.pdb   | 3aligned.xtc
    '''

    kwargs["func"] = lambda args: print(syntax, file=args.fout)

    help = "Print an example of the input file for the align and fit commands."
    parser = subparsers.add_parser("input_syntax", help=help, description=help)

    msg = 'File to write the example syntax to. By default standard output.'
    parser.add_argument('-o', metavar="FNAME", dest='fout',
                        type=argparse.FileType("w"),
                        default=sys.stdout, help=msg)

    parser.set_defaults(**kwargs)

    return parser


def create_align_subparser(subparsers: argparse._SubParsersAction, **kwargs) -> argparse.ArgumentParser:
    """
    Create the subparser for the align command. **kwargs are added to set_defaults
    """
    help = "Align sequences and find common atoms.\nIn other words, find out how to compare apples to oranges."
    parser = subparsers.add_parser("align", help=help, description=help)

    msg = "Input file. For the example syntax run 'apples2apples input_syntax'."
    parser.add_argument('-i', metavar="FNAME", dest='input_file',
                        type=inout.wrap_input_file_type("align"),
                        help=msg, required=True)

    msg = 'Temporary directory for aligment files.'
    parser.add_argument('-t', '--temp', metavar="DIR", dest='temp_directory',
                        type=inout.temporary_directory, help=msg, default=".")

    msg = "Selection for atoms when aligned residues are not the same. " \
          "Options are either the whole backbone \'backbone\' or the " \
          "alpha-carbon atoms \'CA\'. Default is \'CA\'."
    parser.add_argument('-s', metavar="CA|backbone", dest='not_aligned_sel',
                        type=inout.not_aligned_selection,
                        default='ca', help=msg)

    parser.add_argument("-p", "--print-selections",
                        dest="print_sels", action="store_true",
                        help="Print the selections to stdout.")

    parser.set_defaults(**kwargs)

    return parser


def create_fit_subparser(subparsers: argparse._SubParsersAction, **kwargs) -> argparse.ArgumentParser:
    """
    Create the subparser for the fit command. **kwargs are added to set_defaults
    """
    help = "Fit trajectories (in space).\nIn other words, put apples on oranges in preparation for comparison of apples to oranges."
    parser = subparsers.add_parser("fit", help=help, description=help)

    msg = "Input file. For the example syntax run 'apples2apples input_syntax'. " \
          "NOTE that in this case the index file inputs will be read from the 'output_ndx' " \
          "columns, as these are where the align command wrote its output."
    parser.add_argument('-i', metavar="FNAME", dest='input_file',
                        type=inout.wrap_input_file_type("fit"),
                        help=msg, required=True)

    msg = 'Reference frame used for the alignment. The reference is taken from the first pdb listed in the input.'
    parser.add_argument('-r', '--ref', metavar="N", dest='ref_frame',
                        type=inout.temporary_directory, help=msg, default=".")

    parser.set_defaults(ref_frame=0, **kwargs)

    return parser
