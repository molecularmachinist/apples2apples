from typing import List
import numpy as np
import Bio.SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.SeqRecord import SeqRecord

from . import utils


def read_seqs_from_aligment_file(aligned_fasta_file: str, ids: List[str]) -> List[str]:

    seqs = ['' for _ in ids]
    i_seq = -1
    with open(aligned_fasta_file, 'r') as file:
        for line in file:
            if line[:4] == '>pdb':
                i_seq += 1
            else:
                # [:-1] to remove \n from the end of the line
                seqs[i_seq] += line[:-1]
    return seqs


def aligned_sequences(ids: List[str], records: List[SeqRecord], temp: str):

    unaligned_fasta_file = '{}/unaligned.fasta'.format(temp)
    Bio.SeqIO.write(records, unaligned_fasta_file, 'fasta')

    aligned_fasta_file = '{}/aligned.fasta'.format(temp)
    clustalomega_cline = ClustalOmegaCommandline(
        infile=unaligned_fasta_file, outfile=aligned_fasta_file, verbose=True, auto=True, force=True)
    clustalomega_cline()

    return read_seqs_from_aligment_file(aligned_fasta_file, ids)


def find_apples2apples(seqs: List[str], resids: List[List[int]], not_aligned_sel: str) -> List[str]:

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
