import argparse
import os
import pathlib
from typing import Callable, Dict, List, Tuple, Union
import MDAnalysis as mda
from Bio.SeqRecord import SeqRecord


def read_ndx(ndx: str, verbose=False) -> Dict[str, List[int]]:
    # Return empty dictionary if ndx is None
    if (ndx is None):
        return {}
    if verbose:
        print("\nReading index groups from %s" % ndx)
    indexes = {}
    groups = []
    current = None

    with open(ndx) as f:
        for line in f:
            line = line.strip()
            if (line.startswith("[") and line.endswith("]")):
                current = line[1:-1].strip()
                indexes[current] = []
                groups.append(current)
                continue
            parts = line.split()
            for i in parts:
                indexes[current].append(int(i))

    if (verbose):
        print("Found %d groups" % len(groups))
        print("    %-20s%10s\n" % ("Group name", "atoms"))
        for i, g in enumerate(groups):
            print("%2d. %-20s%10d" % (i+1, g, len(indexes[g])))

        print()
    return indexes

###############################################################################################################
# Functions to parse the input file
###############################################################################################################


def is_file_readable(file: str):
    # Try whether the files can be read or not

    try:
        with open(file, 'r'):
            pass
    except Exception as err:
        raise argparse.ArgumentTypeError(err)


def is_file_writable(file: str):
    # Try whether the files can be writen or not

    try:
        with open(file, 'w'):
            pass
    except Exception as err:
        raise argparse.ArgumentTypeError(err)


def checkinput_pdbs_type(arg: str) -> List[str]:

    input_pdbs = arg.split()

    if len(input_pdbs) < 2:

        msg = 'At least two pdb files are required as an input. You provided only {}: {}'.format(
            len(input_pdbs), arg)
        raise argparse.ArgumentTypeError(msg)

    else:
        # Try whether the files can be read or not
        for input_pdb in input_pdbs:
            is_file_readable(input_pdb)

    return input_pdbs


__col_info = {
    "required_align": ['input_pdb', 'output_ndx'],
    "required_fit": ['input_pdb', "input_xtc", 'output_ndx', 'output_traj'],
    "allowed_columns": {'input_pdb', "input_xtc", 'selection',
                        'first_residue_index', 'output_ndx', 'output_traj', 'output_pdb'}
}


def check_required_column_titles_exits(column_titles: str, input_file: str, required: str) -> Tuple[Dict[str, int], int]:

    column_titles = {column_title.strip(): i
                     for i, column_title in enumerate(column_titles.split('|'))}

    n_columns = len(column_titles)
    if n_columns > 6:
        msg = 'Input file \'{}\' is containing {} columns. Input file should contain at most 6 columns'.format(
            input_file, n_columns)
        raise argparse.ArgumentTypeError(msg)

    for required_column_title in __col_info[required]:
        if required_column_title not in column_titles:

            msg = 'Given input file \'{}\' is missing \'{}\' column.'.format(
                input_file, required_column_title)
            raise argparse.ArgumentTypeError(msg)

    if ('output_ndx' not in column_titles) and ('output_pdb' not in column_titles):
        msg = 'Given input file \'{}\' is missing both columns \'output_ndx\' and \'output_pdb\'. At least one of them is required.'.format(
            input_file)
        raise argparse.ArgumentTypeError(msg)

    for title in column_titles:
        if (title not in __col_info["allowed_columns"]):
            msg = f'Given input file \'{input_file}\' has an unrecognized column title \'{title}\'. Allowed values are {__col_info["allowed_columns"]}.'
            raise argparse.ArgumentTypeError(msg)

    return column_titles, n_columns


def split_inputfile(lines: List[str], n_columns: int) -> Tuple[int, List[List[str]]]:
    n_pdbs = len(lines)-1
    if n_pdbs < 2:
        msg = 'At least two pdb files are required as an input. You provided only {}.'.format(
            n_pdbs)
        raise argparse.ArgumentTypeError(msg)

    input_matrix = [[] for _ in range(n_columns)]

    for i, line in enumerate(lines[1:]):
        line_splitted = [columns.strip() for columns in line.split('|')]

        if len(line_splitted) != n_columns:
            msg = 'Given input file contains {} columns, but line of {}th pdb in the given input file contains {} columns.'.format(
                n_columns, i+1, len(line_splitted))
            raise argparse.ArgumentTypeError(msg)

        for j, element in enumerate(line_splitted):
            input_matrix[j].append(element)

    return n_pdbs, input_matrix


def read_columns_and_rows(lines: List[str], input_file: str, funcname: str):

    column_dict, n_columns = check_required_column_titles_exits(
        lines[0], input_file, f"required_{funcname}")

    n_pdbs, input_matrix = split_inputfile(lines, n_columns)

    # check that input files are readable
    input_pdbs = input_matrix[column_dict['input_pdb']]
    for input_pdb in input_pdbs:
        is_file_readable(input_pdb)

    inputdata = {}
    for key in __col_info['allowed_columns']:
        if (key in column_dict):
            inputdata[key] = input_matrix[column_dict[key]]
        else:
            inputdata[key] = None

    inputdata["id"] = []
    for i in range(n_pdbs):
        inputdata["id"].append('pdb{}'.format(i+1))

    # These are needed as dicts
    for key in ("input_xtc", "output_traj"):
        if key in column_dict:
            inputdata[key] = {pdb: trj for pdb, trj in zip(
                input_pdbs, inputdata[key])}

    # set selections to "all" if needed
    if inputdata["selection"] is None:
        inputdata["selection"] = ["all" for p in input_pdbs]

    # create universes and use selection strings
    univs = {}
    records = []
    subsets = []
    for i, (pdb, sel) in enumerate(zip(input_pdbs, inputdata["selection"])):
        if (pdb not in univs):
            univs[pdb] = mda.Universe(pdb)

        u = univs[pdb]

        subset = u.select_atoms(sel)
        if (len(subset) == 0):
            msg = 'Selection of line {} ("{}") is empty'.format(
                i+1, sel)
            raise argparse.ArgumentTypeError(msg)
        subsets.append(subset)

        record = subset.residues.sequence(id=inputdata["id"][i])
        records.append(record)

    inputdata["univ"] = univs
    inputdata["subset"] = subsets
    inputdata["record"] = records

    # convert first_residue_indexes to integers or get from selections
    if ("first_residue_index" in column_dict):
        for i in range(n_pdbs):
            try:
                inputdata["first_residue_index"][i] = int(
                    inputdata["first_residue_index"][i])
            except Exception as err:
                msg = 'Error occurred in \'first_residue_index\' of line {}: {}'.format(
                    i+1, err)
                raise argparse.ArgumentTypeError(msg)
    else:
        inputdata["first_residue_index"] = [sel[0].resid for sel in subsets]

    return inputdata


def input_file_type(input_file: str, command: str) -> dict:

    is_file_readable(input_file)

    with open(input_file) as file:

        lines = []
        for line in file:
            if line.strip() and (not line.lstrip().startswith("#")):
                lines.append(line)

    return read_columns_and_rows(lines, input_file, command)


def wrap_input_file_type(command: str) -> Callable[[str], dict]:
    return lambda input_file: input_file_type(input_file, command)


def temporary_directory(temp: str):
    temp_path = pathlib.Path(temp)
    try:
        temp_path.mkdir(exist_ok=True)
    except Exception as err:
        msg = 'Error occurred while trying to make sure temporary directory exists: {}\nMake sure its parents exist and are writable'.format(
            err)
        raise argparse.ArgumentTypeError(msg)

    test_file = temp_path / f'apples2apples_test_file_{os.getpid()}.txt'
    try:
        with test_file.open('w'):
            pass
        os.remove(test_file)
    except Exception as err:
        msg = 'Error occurred while trying to write and remove a test file in the given temporary directory: {}'.format(
            err)
        raise argparse.ArgumentTypeError(msg)

    return temp


def not_aligned_selection(sel: str):
    if (sel != 'backbone') and (sel != 'ca') and (sel != 'CA'):
        msg = 'Options are \'backbone\',  \'CA\' and \'ca\'. You provided: {}'.format(
            sel)
        raise argparse.ArgumentTypeError(msg)

    if sel == 'backbone':
        return 'backbone'
    else:
        return 'name CA'
