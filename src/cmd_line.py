import argparse
from typing import Tuple

from . import inout


def create_main_parser() -> Tuple[argparse.ArgumentParser, argparse._SubParsersAction]:
    description = "apples2apples: A tool for comparing apples to oranges."
    epilog = "For more detailed explanations, see the README.md"
    parser = argparse.ArgumentParser(
        prog="apples2apples", description=description, epilog=epilog)

    msg = "Specify what the program should do. Only one command should be given per run. 'apples2apples <cmd> -h' to get command specific option."
    subparsers = parser.add_subparsers(description=msg)
    subparsers.required = True
    subparsers.dest = 'command'

    return parser, subparsers


def create_input_syntax_subparser(subparsers: argparse._SubParsersAction, **kwargs) -> argparse.ArgumentParser:
    msg = '''
# Input file example.
# The symbol # is used for comment lines.
# The first non-empty line after these comments includes column titles.
# The columns are separated from eachother by pipe symbol |.
# The columns input_pdbs, selections, first_residue_index and either output_ndx or output_pdb
# are required. In the input_pdbs at least 2 pdb files are required.
# See  https://docs.mdanalysis.org/stable/documentation_pages/selections.html
# for syntax of the selections column. The first_residue_index column contains
# the residue index used in the pdb file of corresponding line.
# The columns output_ndx and output_pdb are to specify output files
# Note that the pipe symbols don't have to be aligned.

input_pdbs |  selections              | first_residue_index | output_ndx | output_pdb
1.pdb      |  segid A and resid 40:60 | 40                  | 1.ndx      | 1out.pdb
2.pdb      |  segid B and resid 10:55 | 10                  | 2.ndx      | 2out.pdb
3.pdb      |  protein                 | 483                 | 3.ndx      | 3out.pdb
    '''

    kwargs["func"] = lambda args: print(msg)

    help = "Print an example of the input file for the align and model commands."
    parser = subparsers.add_parser("input_syntax", help=help, description=help)

    parser.set_defaults(**kwargs)

    return parser


def create_align_subparser(subparsers: argparse._SubParsersAction, **kwargs) -> argparse.ArgumentParser:
    help = "Align sequences and find common atoms.\nIn other words, find out how to compare apples to oranges."
    parser = subparsers.add_parser("align", help=help, description=help)

    msg = "Input file. For the example syntax run 'apples2apples input_syntax'."
    parser.add_argument('-i', dest='input_file',
                        type=inout.input_file_type, help=msg, required=True)

    msg = 'temporary directory for aligment files.'
    parser.add_argument('-t', '--temp', dest='temp_directory',
                        type=inout.temporary_directory, help=msg, required=True)

    msg = 'Selection for atoms when aligned residues are not the same. Options are either the whole backbone \'backbone\' or the alpha-carbon atoms \'ca\'. Default is \'ca\'.'
    parser.add_argument('-s', dest='not_aligned_sel',
                        type=inout.not_aligned_selection, default='ca', help=msg)

    parser.set_defaults(**kwargs)

    return parser


def create_model_subparser(subparsers: argparse._SubParsersAction, **kwargs) -> argparse.ArgumentParser:
    help = "Fit sequences (in space) and train ML models.\nIn other words, compare apples to oranges."
    parser = subparsers.add_parser("model", help=help, description=help)

    msg = "Input file. For the example syntax run 'apples2apples input_syntax'. " \
          "NOTE that in this case the real inputs will be read from the 'output' columns, " \
          "as these are where the align command wrote its output."
    parser.add_argument('-i', dest='input_file',
                        type=inout.input_file_type, help=msg, required=True)

    parser.set_defaults(**kwargs)

    return parser
