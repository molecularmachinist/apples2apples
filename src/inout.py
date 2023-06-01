import argparse
import os
import pathlib
from typing import Callable, Dict, List, Tuple, Union
import MDAnalysis as mda
from Bio.SeqRecord import SeqRecord


def read_ndx(ndx: pathlib.Path, verbose=False) -> Dict[str, List[int]]:
    """
    Read the gromacs .ndx file

    Parameters:
    -----------
    ndx: str
        The path to the index file
    verbose: bool=False
        whether to print the found group names and sizes

    Returns:
    --------
    indexes: dict[str, list[int]]
        A dictionary of the index groups as lists of string, with the group names as keys
    """
    # Return empty dictionary if ndx is None
    if (ndx is None):
        return {}
    if verbose:
        print("\nReading index groups from %s" % ndx)
    indexes = {}
    groups = []
    current = None

    with ndx.open() as f:
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


def is_file_readable(file: pathlib.Path):
    """
    Make sure the file is readable. Raises an ArgumentTypeError if not.
    Returns the input parameter unchanged.
    """
    try:
        with file.open('r'):
            pass
    except Exception as err:
        raise argparse.ArgumentTypeError(err)
    return file


def is_file_writable(file: pathlib.Path, ):
    """
    Make sure the file is writable. Raises an ArgumentTypeError if not.
    Returns the input parameter unchanged.
    """
    try:
        file.parent.mkdir(parents=True, exist_ok=True)
        with file.open('w'):
            pass
    except Exception as err:
        raise argparse.ArgumentTypeError(err)
    return file


__col_info = {
    "required_align": ['input_pdb', 'output_ndx'],
    "required_fit": ['input_pdb', "input_xtc", 'output_ndx', 'output_traj'],
    "allowed_columns": {'input_pdb', "input_xtc", 'selection',
                        'output_ndx', 'output_traj', 'output_pdb'},
    "used_align": {'input_pdb', 'selection', 'output_ndx', 'output_pdb'},
    "used_fit": {'input_pdb', "input_xtc", 'output_ndx', 'output_traj'},
    "paths": {"input_pdb", "input_xtc", "output_ndx", "output_pdb", "output_traj"}
}


def check_required_column_titles_exits(column_titles: str, input_file: pathlib.Path, required: str) -> Tuple[Dict[str, int], int]:
    """
    Checks column titles and returns a mapping of title to index.
    Raises an ArgumentTypeError if not required title is missing or an unknown one is found.

    Parameters:
    -----------
    column_titles: str
        The first nonempty and noncomment line of the input file. Should have the column titles.
    input_file: Path
        Path to the input file.
    required: str
        The key to get the required column name field.
        Should be "required_align" or "required_fit" depending on the running command.

    Returns:
    --------
    column_titles: dict[str, int]
        A dictionary mapping the values of column names to column indexes
    n_columns: int
        The number of columns found
    """

    column_titles = {column_title.strip(): i
                     for i, column_title in enumerate(column_titles.split('|'))}

    n_columns = len(column_titles)
    if n_columns > 6:
        msg = f'Input file \'{input_file}\' is containing {n_columns} columns. Input file should contain at most 6 columns'
        raise argparse.ArgumentTypeError(msg)

    for required_column_title in __col_info[required]:
        if required_column_title not in column_titles:

            msg = f'Given input file \'{input_file}\' is missing \'{required_column_title}\' column.'
            raise argparse.ArgumentTypeError(msg)

    if ('output_ndx' not in column_titles) and ('output_pdb' not in column_titles):
        msg = f'Given input file \'{input_file}\' is missing both columns \'output_ndx\' and \'output_pdb\'. At least one of them is required.'
        raise argparse.ArgumentTypeError(msg)

    for title in column_titles:
        if (title not in __col_info["allowed_columns"]):
            msg = f'Given input file \'{input_file}\' has an unrecognized column title \'{title}\'. Allowed values are {__col_info["allowed_columns"]}.'
            raise argparse.ArgumentTypeError(msg)

    return column_titles, n_columns


def split_inputfile(lines: List[str], n_columns: int) -> Tuple[int, List[List[str]]]:
    """
    Split the input file into a "matrix" of list of lists of strings, where the first index is over the columns and the second over rows.

    Parameters:
    -----------
    lines: list[str]
        List of the lines in the input file.
    n_columns: int
        Number of columns in the file.

    Returns:
    --------
    n_pdbs: int
        number of pdb files found (number of rows)
    input_matrix: list[list[str]]
        The input file split into a matrix of strings
    """
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


def read_columns_and_rows(lines: List[str], input_file: pathlib.Path, funcname: str):
    """
    Read the input file content and return a dictionary of the data.

    Parameters:
    -----------
    lines: list[str]
        List of the lines in the input file.
    input_file: Path
        Path to the input file.
    funcname: str
        The name of the running command, should be "align" or "fit".

    Returns:
    --------
    inputdata: dict
        The dictionary with the input data
    """

    column_dict, n_columns = check_required_column_titles_exits(lines[0],
                                                                input_file,
                                                                f"required_{funcname}")

    n_pdbs, input_matrix = split_inputfile(lines, n_columns)

    inputdata = {}
    for key in __col_info[f"used_{funcname}"]:
        if (key in column_dict):
            inputdata[key] = input_matrix[column_dict[key]]
            if (key in __col_info["paths"]):
                inputdata[key] = [pathlib.Path(item)
                                  for item in inputdata[key]]
        else:
            inputdata[key] = None

    # check that input files are readable
    input_pdbs = [is_file_readable(input_pdb)
                  for input_pdb in inputdata["input_pdb"]]

    inputdata["id"] = []
    for i in range(n_pdbs):
        inputdata["id"].append('pdb{}'.format(i+1))

    # These are needed as dicts
    for key in ("input_xtc", "output_traj"):
        if key in inputdata:
            inputdata[key] = {pdb: trj for pdb, trj in
                              zip(input_pdbs, inputdata[key])}

    # set selections to "all" if needed
    if "selection" in inputdata and inputdata["selection"] is None:
        inputdata["selection"] = ["all" for p in input_pdbs]

    # create universes and use selection strings
    univs = {}
    for pdb in input_pdbs:
        if pdb not in univs:
            univs[pdb] = mda.Universe(str(pdb))

    inputdata["univ"] = univs

    if (funcname == "align"):
        records = []
        subsets = []
        for i, (pdb, sel) in enumerate(zip(input_pdbs, inputdata["selection"])):
            u = univs[pdb]

            subset = u.select_atoms(sel)
            if (len(subset) == 0):
                msg = f'Selection of line {i+1} ("{sel}") is empty'
                raise argparse.ArgumentTypeError(msg)
            subsets.append(subset)

            record = subset.residues.sequence(id=inputdata["id"][i])
            records.append(record)

        inputdata["subset"] = subsets
        inputdata["record"] = records

    return inputdata


def input_file_type(input_file: pathlib.Path, command: str) -> dict:
    """
    A file "type" for argparse. Reads and checks the given input file,
    returning its contents as a dictionary.
    """

    is_file_readable(input_file)

    with input_file.open() as file:

        lines = []
        for line in file:
            if line.strip() and (not line.lstrip().startswith("#")):
                lines.append(line)

    return read_columns_and_rows(lines, input_file, command)


def wrap_input_file_type(command: str) -> Callable[[str], dict]:
    """
    A wrapper for the input_file_type, allowing to specify the command.
    E.g. to get the input for fit, give wrap_input_file_type("fit") as a type
    to argparse.
    """
    return lambda input_file: input_file_type(pathlib.Path(input_file), command)


def temporary_directory(temp: str):
    """
    Check that the directory is writable and raise an exception if not.
    Returns the argument if successfull.
    """
    temp_path = pathlib.Path(temp)
    try:
        temp_path.mkdir(exist_ok=True)
    except Exception as err:
        msg = f'Error occurred while trying to make sure temporary directory exists: {err}\n' \
            'Make sure its parents exist and is writable'
        raise argparse.ArgumentTypeError(msg)

    test_file = temp_path / f'apples2apples_test_file_{os.getpid()}.txt'
    try:
        with test_file.open('w'):
            pass
        os.remove(test_file)
    except Exception as err:
        msg = 'Error occurred while trying to write and remove a test file in the given temporary directory: {err}'
        raise argparse.ArgumentTypeError(msg)

    return temp_path


def not_aligned_selection(sel: str):
    """
    Check that the sel is one of the two allowed possibilities, ("CA" or "backbone").
    """
    if (sel != 'backbone') and (sel != 'ca') and (sel != 'CA'):
        msg = 'Options are \'backbone\',  \'CA\' and \'ca\'. You provided: {}'.format(
            sel)
        raise argparse.ArgumentTypeError(msg)

    if sel == 'backbone':
        return 'backbone'
    else:
        return 'name CA'
