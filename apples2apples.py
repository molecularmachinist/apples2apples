import argparse
import MDAnalysis as mda
import numpy as np
import Bio.SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from warnings import catch_warnings, simplefilter





###############################################################################################################
#################################  Functions to parse the input file
###############################################################################################################

def is_file_readable(file):
    # Try whether the files can be read or not

    try:
        with open(file, 'r'):
            pass
    except Exception as err:
        raise argparse.ArgumentTypeError(err)


def is_file_writable(file):
    # Try whether the files can be writen or not

    try:
        with open(file, 'w'):
            pass
    except Exception as err:
        raise argparse.ArgumentTypeError(err)



def checkinput_pdbs_type(arg):
    
    input_pdbs  =arg.split()

    if len(input_pdbs) < 2:

        msg = 'At least two pdb files are required as an input. You provided only {}: {}'.format(len(input_pdbs),arg)
        raise argparse.ArgumentTypeError(msg)

    else:
        # Try whether the files can be read or not
        for input_pdb in input_pdbs:
            is_file_readable(input_pdb)

    return input_pdbs




def check_required_column_titles_exits(column_titles, input_file):


    column_titles_splitted = [column_title.strip() for column_title in column_titles.split('|')]
    
    n_columns = len(column_titles_splitted)
    if n_columns > 5:
        msg = 'Input file \'{}\' is containing {} columns. Input file should contain at most 5 columns'.format(input_file, n_columns)
        raise argparse.ArgumentTypeError(msg)


    for required_column_title in ['input_pdbs', 'selections','first_residue_index']:
        if required_column_title not in column_titles_splitted:

            msg = 'Given input file \'{}\' is missing \'{}\' column.'.format(input_file, required_column_title)
            raise argparse.ArgumentTypeError(msg)

    if ('output_ndx' not in  column_titles_splitted) and ('output_pdb' not in column_titles_splitted):
        msg = 'Given input file \'{}\' is missing both columns \'output_ndx\' and \'output_pdb\'. At least one of them is required.'.format(input_file)
        raise argparse.ArgumentTypeError(msg)


    return dict( zip( range(n_columns), column_titles_splitted  )  ), n_columns
    



def read_columns_and_rows(lines, input_file):

    column_dict, n_columns = check_required_column_titles_exits(lines[0], input_file)

    column_dict_inv = {v: k for k, v in column_dict.items()}

    n_pdbs = len(lines)-1
    if n_pdbs < 2:
        msg = 'At least two pdb files are required as an input. You provided only {}.'.format(n_pdbs)
        raise argparse.ArgumentTypeError(msg)


    input_matrix = [[] for _ in range(n_columns)]




    for i, line in enumerate(lines[1:]):
        line_splitted = [columns.strip() for columns in line.split('|')]

        if len(line_splitted) != n_columns:
            msg = 'Given input file \'{}\' contains {} columns, but line of {}th pdb in the given input file contains {} columns.'.format(input_file, n_columns, i+1, len(line_splitted))
            raise argparse.ArgumentTypeError(msg)

        for j, element in enumerate(line_splitted):
            input_matrix[j].append(element)


    # check that input files are readable
    input_pdbs = input_matrix[  column_dict_inv['input_pdbs']]
    for input_pdb in  input_pdbs:
        is_file_readable(input_pdb)

    # check that output files are readable
    if 'output_ndx' in column_dict_inv.keys():
        output_ndxs = input_matrix[  column_dict_inv['output_ndx']]
        for ndx in output_ndxs:
            is_file_writable(ndx)

    else:
        output_ndxs = None

    if 'output_pdb' in column_dict_inv.keys():
        output_pdbs  = input_matrix[  column_dict_inv['output_pdb']]
        for pdb in output_pdbs:
            is_file_writable(pdb)
    else:
        output_pdbs = None

    # convert first_residue_indexes to integers
    first_residue_indexes = input_matrix[  column_dict_inv['first_residue_index']]
    for i in range(n_pdbs):
        try:
            
            input_matrix[column_dict_inv['first_residue_index']][i] = int(first_residue_indexes[i])
        except Exception as err:
            msg = 'Error occurred in \'first_residue_index\' of line {}: {}'.format(i+1,err)
            raise argparse.ArgumentTypeError(msg)


    # creat universes and use selection strings

    ids = []
    for i in range(n_pdbs):
        ids.append('pdb{}'.format(i+1))

    universes = []
    records = []
    subsets = []

    i = 0
    selections = input_matrix[  column_dict_inv['selections']]
    for pdb, sel in zip(input_pdbs,selections):
        u = mda.Universe(pdb)
        
        universes.append(u)
        
        subset = u.select_atoms(sel)
        subsets.append(subset)
        
        record = subset.residues.sequence(id=ids[i])
        records.append(record)
        
        i += 1



    return ids, universes, records, subsets, first_residue_indexes, output_ndxs, output_pdbs
        



def input_file_type(input_file):

    is_file_readable(input_file)

    with open(input_file) as file:

        lines = []
        for line in file:
            if line.strip() and (line.lstrip()[0] != "#"):
                lines.append(line)


    return read_columns_and_rows(lines, input_file)


def temporary_directory(temp):

    test_file = '{}/test_file.txt'.format(temp)
    try:
        with open(test_file, 'w'):
            pass
    except Exception as err:
        msg = 'Error occurred while trying to write a test file to the given temporary directory: {}'.format(err)
        raise argparse.ArgumentTypeError(msg)

    return temp

def not_aligned_selection(sel):
    if (sel != 'backbone') and (sel != 'ca'):
        msg = 'Options are \'backbone\' and \'ca\'. You provided: {}'.format(sel)
        raise argparse.ArgumentTypeError(msg)

    if sel == 'ca':
        return 'name CA'
    else:
        return 'backbone'


def create_parser():


    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    # input pdb files
    #msg = 'String of two or more pdb files as an input, e.g. -i "1.pdb 2.pdb 3.pdb"'
    msg = '''Input file with following syntax:

    # Input file example.
    # The symbol # is used for comment lines.
    # The first non-empty line after these comments includes column titles.
    # The columns are separated from eats other by pipe symbol |.
    # The columns input_pdbs, selections, first_residue_index and either output_ndx or output_pdb
    # are required. In the input_pdbs at least 2 pdb files are required.
    # See  https://docs.mdanalysis.org/stable/documentation_pages/selections.html
    # for syntax of the selections column. The first_residue_index column contains
    # the residue index used in the pdb file of corresponding line.
    # The columns output_ndx and output_pdb are to specify output files
    # Note that the pipe symbols don't have to be align.

    input_pdbs |  selections              | first_residue_index | output_ndx | output_pdb
    1.pdb      |  segid A and resid 40:60 | 40                  | 1.ndx      | 1out.pdb
    2.pdb     |  segid B and resid 10:55 | 10                  | 2.ndx      | 2out.pdb
    3.pdb      |  protein                 | 483                 | 3.ndx      | 3out.pdb
    '''
    parser.add_argument('-i', dest='input_file', type=input_file_type, help=msg,required=True)


    msg = 'temporary directory for aligment files.'
    parser.add_argument('-t','--temp', dest='temp_directory', type=temporary_directory, help=msg,required=True)

    
    msg = 'Selection for atoms when aligned residues are not the same. Options are either the whole backbone \'backbone\' or the alpha-carbon atoms \'ca\'. Default is \'ca\'.' 
    parser.add_argument('-s', dest='not_aligned_sel', type=not_aligned_selection, default='ca', help=msg)

    return parser


###############################################################################################################
#################################  To do the aligment and apples2apples
###############################################################################################################

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



def main():
    parser  = create_parser()

    args = parser.parse_args()



    ids, universes, records, subsets, first_residue_indexes, output_ndxs, output_pdbs = args.input_file 
    temp = args.temp_directory
    not_aligned_sel = args.not_aligned_sel


    seqs = aligned_sequences(ids, records, temp)


    indeces = [i-1 for i in first_residue_indexes]

    sels = find_apples2apples(seqs, indeces, not_aligned_sel)



    common_sels = []
    for sel,subset in zip(sels,subsets):
        common_sel = subset.select_atoms(sel)
        common_sels.append(common_sel)

    if output_ndxs is not None:
        for common_sel, ndx_file in zip(common_sels, output_ndxs):
            with mda.selections.gromacs.SelectionWriter(ndx_file, mode='w') as ndx:
                ndx.write(common_sel, name='apples2apples')

    if output_pdbs is not None:
        with catch_warnings():
            simplefilter("ignore", category=(DeprecationWarning,UserWarning))

            for common_sel, pdb in zip(common_sels, output_pdbs):
                common_sel.write(pdb)



if __name__ == '__main__':


    main()
