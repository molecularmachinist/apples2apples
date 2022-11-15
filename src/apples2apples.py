import pathlib
from typing import List
import numpy as np
import Bio.SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.SeqRecord import SeqRecord

from . import utils


def read_seqs_from_aligment_file(aligned_fasta_file: pathlib.Path) -> List[str]:
    """
    Reads sequences from 'aligned_fasta_file'. Returns a list of each sequence as a single string.
    """
    seqs = []
    with aligned_fasta_file.open('r') as file:
        for line in file:
            if line[:4] == '>pdb':
                seqs.append("")
            else:
                # strip to remove \n from the end of the line
                seqs[-1] += line.strip()
    return seqs


def aligned_sequences(records: List[SeqRecord], temp: pathlib.Path):
    """
    Takes the list of the sequences, aligns them and returns the alignment results as a list of strings.

    Parameters
    ----------
    records: list[SeqRecord]
        A list of SeqRecord objects, which are the sequences to align
    temp: str
        The path to a directory to use for the temporary fasta files.
        The fasta files will be named unaligned.fasta and aligned.fast
        and will not be removed when done.

    Returns
    -------
    sequences: list[str]
        list of strings, where each string is a single sequence after the alignment.
    """

    unaligned_fasta_file = temp / "unaligned.fasta"
    Bio.SeqIO.write(records, unaligned_fasta_file, 'fasta')

    aligned_fasta_file = temp / "aligned.fasta"
    clustalomega_cline = ClustalOmegaCommandline(infile=unaligned_fasta_file,
                                                 outfile=aligned_fasta_file,
                                                 verbose=True,
                                                 auto=True,
                                                 force=True)
    clustalomega_cline()

    return read_seqs_from_aligment_file(aligned_fasta_file)


def find_apples2apples(seqs: List[str], resids: List[List[int]], not_aligned_sel: str) -> List[str]:
    """
    Find the MDAnalysis selection strings from the aligned sequences, which will results in atom
    selections of the same size.

    Parameters:
    -----------
    seqs: list[str]
        List of the aligned sequences, each as a single string
    resids: list[list[int]]
        List of the lists of resids for each sequence
    not_aligned_sel: str
        A string to use as a selection for the residues that are aligned, but do not match

    Returns:
    --------
    sels: list[str]
        list of the selection strings for each sequence
    """

    indeces = np.zeros(len(seqs), dtype=int)

    backbones = [[] for _ in seqs]
    whole_ress = [[] for _ in seqs]

    for aas in zip(*seqs):
        aas = np.asarray(aas)

        bools = aas != '-'

        if np.all(bools):

            if np.all(aas == aas[0]):
                for i, index in enumerate(indeces):
                    whole_ress[i].append(resids[i][index])

            else:
                for i, index in enumerate(indeces):
                    backbones[i].append(resids[i][index])

            indeces = indeces + 1
        else:
            indeces += bools

    print(f"Found {len(whole_ress[0])} common residues "
          f"and {len(backbones[0])} alignable but different")
    sels = []
    for backbone, whole_res in zip(backbones, whole_ress):
        if not backbone:
            if not whole_res:
                sel = ''
            else:
                sel = 'resid {}'.format(
                    utils.compressed_resid_list(whole_res))

        elif not whole_res:
            sel = '{} and resid {}'.format(
                not_aligned_sel,
                utils.compressed_resid_list(backbone))
        else:
            sel = '({} and resid {}) or (resid {})'.format(
                not_aligned_sel,
                utils.compressed_resid_list(backbone),
                utils.compressed_resid_list(whole_res))

        sels.append(sel)

    return sels
