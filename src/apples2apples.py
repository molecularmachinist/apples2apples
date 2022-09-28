#!/usr/bin/env python3
import numpy as np
import Bio.SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline


def read_seqs_from_aligment_file(aligned_fasta_file, ids):
    
    
    seqs = ['' for _ in ids]

    
    
    with open(aligned_fasta_file,'r') as file:
        for line in file:
            if line[:4]== '>pdb':
                i_seq = int(line[4])-1
            else:
                seqs[i_seq] += line[:-1] # [:-1] to remove \n from the end of the line
    return seqs



def aligned_sequences(ids, records, temp):

    unaligned_fasta_file = '{}/unaligned.fasta'.format(temp)
    Bio.SeqIO.write(records, unaligned_fasta_file, 'fasta')


    aligned_fasta_file = '{}/aligned.fasta'.format(temp)
    clustalomega_cline = ClustalOmegaCommandline(infile=unaligned_fasta_file, outfile=aligned_fasta_file, verbose=True, auto=True, force=True)
    clustalomega_cline()

    return read_seqs_from_aligment_file(aligned_fasta_file, ids)





def find_apples2apples(seqs, indeces, not_aligned_sel):
    
    indeces = np.asarray(indeces)
    
    backbones = ['' for _ in seqs]
    whole_ress = ['' for _ in seqs]
    

    for aas in zip(*seqs):
        aas = np.asarray(aas)

        
        bools = np.array([aa != '-' for aa in aas])
        

        
        if np.all(bools):
            

            indeces = indeces + 1
            
            if np.all( aas==aas[0] ):
                for i, index in enumerate(indeces):
                    whole_ress[i] += ' {}'.format(index)

            else:
                for i, index in enumerate(indeces):
                    backbones[i] += ' {}'.format(index)
                
        
        else:
            for i,booL in enumerate(bools):
                if booL:
                    indeces[i] += 1

    sels = []
    for backbone, whole_res in zip(backbones, whole_ress):
        if backbone == '':
            if whole_res == '':
                sel = ''
            else:
                sel = 'resid {}'.format(whole_res)

        elif whole_res == '':
            sel = '{} and resid {}'.format(not_aligned_sel, backbone)
        else:
            sel = '({} and resid {}) or (resid {})'.format(not_aligned_sel, backbone,whole_res)
            #sel = 'resid {}'.format(whole_res)
            
        sels.append(sel)


    return sels
