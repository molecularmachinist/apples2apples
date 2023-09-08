import argparse
import sys
from typing import Tuple

from . import inout


def create_main_parser() -> Tuple[argparse.ArgumentParser,
                                  "argparse._SubParsersAction[argparse.ArgumentParser]"]:
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


def create_align_subparser(subparsers: "argparse._SubParsersAction[argparse.ArgumentParser]",
                           **kwargs) -> argparse.ArgumentParser:
    """
    Create the subparser for the align command. **kwargs are added to set_defaults
    """
    help = "Align sequences and find common atoms.\nIn other words, find out how to compare apples to oranges."
    parser = subparsers.add_parser("align", help=help, description=help)

    msg = "Input pdb files, one for each sequence. Specify duplicates to read from same file (remember to set selection)."
    parser.add_argument('-i', "--input-pdbs", metavar="*.pdb", dest='input_pdbs',
                        type=inout.input_pdb_type, nargs="+",
                        help=msg, required=True)

    msg = "Selections, one for each sequence Or a single one to use for all pdbs."
    parser.add_argument('-s', "--selections", metavar="sel", dest='selections',
                        type=str, nargs="*",
                        help=msg, default=["all"])

    msg = "Output indexes. Must have same number of files as input pdbs, or be left out. By default will be numbered 1.ndx, 2.ndx, etc."
    parser.add_argument('-o', "--ouput-indexes", metavar="*.ndx", dest='ouput_ndxs',
                        type=inout.output_ndx_type, nargs="*",
                        help=msg, default=None)

    msg = 'Temporary directory for aligment files.'
    parser.add_argument('-t', '--temp', metavar="DIR", dest='temp_directory',
                        type=inout.temporary_directory, help=msg, default=".")

    msg = "Selection for atoms when aligned residues are not the same. " \
          "Options are either the whole backbone \'backbone\' or the " \
          "alpha-carbon atoms \'CA\'. Default is \'CA\'."
    parser.add_argument('-n', "--not-aligned-sel", metavar="CA|backbone", dest='not_aligned_sel',
                        type=inout.not_aligned_selection,
                        default='ca', help=msg)

    parser.add_argument("-p", "--print-selections",
                        dest="print_sels", action="store_true",
                        help="Print the selections to stdout.")

    msg = "Use \"resid\" instead of \"resindex\" as the selection token. " \
          "This results in a selection string, whcih is more portable to other systems" \
          " as it uses the residue numbers instead of internal indices. Do not use when " \
          "selections include duplicate residue numbers."
    parser.add_argument('--resid', action="store_true", help=msg)

    msg = "Do not make the parent directories for output files even if they don't exist. " \
        "Does not affect the temporary directory, which will always be made."
    parser.add_argument("--no-mkdir",
                        dest="make_dirs", action="store_false",
                        help=msg)

    parser.set_defaults(check_func=check_align_options, **kwargs)

    return parser


def check_align_options(args: argparse.Namespace):
    if (len(args.input_pdbs) < 2):
        raise argparse.ArgumentTypeError(
            "At least two sequences should be specified, but only one was given."
        )
    if (len(args.selections) > 1 and len(args.selections) != len(args.input_pdbs)):
        raise argparse.ArgumentTypeError(
            f"Number of selections ({len(args.selections)}) does not match number of pdbs ({len(args.input_pdbs)})"
        )
    if ((args.ouput_ndxs is not None) and len(args.ouput_ndxs) != len(args.input_pdbs)):
        raise argparse.ArgumentTypeError(
            f"Number of output indexes ({len(args.ouput_ndxs)}) does not match number of pdbs ({len(args.input_pdbs)})"
        )


def create_fit_subparser(subparsers: "argparse._SubParsersAction[argparse.ArgumentParser]",
                         **kwargs) -> argparse.ArgumentParser:
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

    parser.add_argument("--no-mkdir",
                        dest="make_dirs", action="store_false",
                        help="Do not make the parent directories for output files even if they don't exist.")

    parser.set_defaults(ref_frame=0, **kwargs)

    return parser
