# apples2apples

A tool for comparing apples to oranges.

## Requirements

1. Python 3.7 or higher (>=3.10 recommended)
1. NumPy
1. [MDAnalysis](https://docs.mdanalysis.org/stable/index.html)
1. [Biopython](https://biopython.org/) (already a requirement for MDAnalysis)



Either install them manually, or follow the instructions below to let pip take care of that for you.

## Installation


Just run
```sh
pip install .
```
in the project root. This will install all dependencies and install the package. The command line tool should now be runnable as 
```sh
apples2apples [command] [options]
```

Alternatively, you can make manually sure the requirements are installed, and then just run the tool without instalilng as
```sh
/path/to/project/apples2apples/apples2apples [command] [options]
```

## Usage

### Alignment

To start off, you'll need (at least) two structure files for the sequences to compare. If you want to make trajectories of the matching atomsgouprs, fitted for optimal tranlation and rotation, you should of course also have the original trajectories.

Next, you'll need to write an input file, which tells the program what to do. To get an example of the input syntax, you can run 
```sh
apples2apples input_syntax
```
which will print a basic example to stdout.

Let's take a look at a simple two-chain comparison, of the systems included in `input_files`. Here is the contents of the example input included in the project root:

```
 input_pdb            | input_xtc            | selection | output_ndx            | output_traj
 input_files/ampa.pdb | input_files/ampa.xtc | chainID A | output_files/ampa.ndx | output_files/ampa_fitted.xtc
 input_files/nmda.pdb | input_files/nmda.xtc | chainID B | output_files/nmda.ndx | output_files/nmda_fitted.xtc
```
All lines starting with a `#` are considered comments and those line, as well as empty lines are ignored. After these, the first line with content is considered to be the column titles. Column order is not important, as each column is identified by its name. Columns are separated by the pip sign. For the sequence alignment, the `input_pdb` and `output_ndx` are the only required fields. If the selection is not give, it is assumed to be `"all"`. The `output_traj` is not used for the first step, only for the spatial fitting. Another optional field that could be added is the `output_pdb`, where the resulting structure of the selection after alignments would be written. The `input_xtc` is not used for this command yet.

The tool can be run as 
```sh
apples2apples align -i infile.txt -t Temp
```
where the compulsory command line option `-i` specifies the input file and the optional option `-t` the "temporary" directory, where the fasta files wil be written. If the temporary directory is not given, they will be written in the current directory.

For each line after the column titles, the program reads in the structure from `input_pdb` and makes a selection using the selection string in the `selection` column. These should be valid [MDAnalysis selection strings](https://docs.mdanalysis.org/stable/documentation_pages/selections.html). The selections are converted to aminoacid sequences and written to `unaligned.fasta`. Alignment is run on the fasta file, and the result is stored in `aligned.fasta`. After this the program reads the aligned sequences and constructs selections, where matching residues are chosen as a whole, and nonmatching ones either just the alpha carbons or the backbone. The selections are written to `output_ndx` and optionally to `output_pdb`.


### Fitting

The next logical step is to extract the selections form the original tarjectories and fit them for optimal RMSD. The same input file as before can be reused, where also the `input_xtc` column is now compulsory. The `selection` and `output_pdb` are ingored. Note, that here the `output_ndx` is actually used as input, but is named as such, to allow reusing the input file from the previous command.

The next command can be run as
```sh
apples2apples fit -i infile.txt
```

The program reads in the `input_pdb` and selections in `output_ndx` and, using the first frame of the first structure as the reference, translates and rotates the structures to minimise RMSD. The trajectories in `input_xtc` are then read and rewritten to `output_traj`. Also a pdb and index file is made for each new trajectory, with the suffix changed to `.ndx` or `.pdb` (the reason for this is outlined below in the advanced usage part).


### Advanced usage

Reading in one chain per system is as simple as seen above, even for more than two systems. A bit more copmlex situation arises when we want to extract multiple chains from single trajectories. For example we can the the AMPA and NMDA receptors, which both have four chains. We want to make the sequence alignment separately for each chain, but wan the final fitting and the structures as single files. This is why the `fit` command considers any chains coming from the same pdb file as being from the same system. The selection that is fitted and written to `output_traj` is made by concatenating the selection from the different chains in the order they are in the input file.

For example let's consider this input file
```
 input_pdb            | input_xtc            | selection | output_ndx              | output_traj

 input_files/ampa.pdb | input_files/ampa.xtc | chainID A | output_files/ampa_A.ndx | output_files/ampa_fitted.xtc
 input_files/ampa.pdb | input_files/ampa.xtc | chainID B | output_files/ampa_B.ndx | output_files/ampa_fitted.xtc
 input_files/ampa.pdb | input_files/ampa.xtc | chainID C | output_files/ampa_C.ndx | output_files/ampa_fitted.xtc
 input_files/ampa.pdb | input_files/ampa.xtc | chainID D | output_files/ampa_D.ndx | output_files/ampa_fitted.xtc

 input_files/nmda.pdb | input_files/nmda.xtc | chainID B | output_files/nmda_B.ndx | output_files/nmda_fitted.xtc
 input_files/nmda.pdb | input_files/nmda.xtc | chainID C | output_files/nmda_C.ndx | output_files/nmda_fitted.xtc
 input_files/nmda.pdb | input_files/nmda.xtc | chainID D | output_files/nmda_D.ndx | output_files/nmda_fitted.xtc
 input_files/nmda.pdb | input_files/nmda.xtc | chainID A | output_files/nmda_A.ndx | output_files/nmda_fitted.xtc
```

The alignment will be made separately on each chain, so that all chains will share the same number of atoms. However in the fitting command all four AMPA receptor chains are considered as a part of a single system (meaning that only the last `input_xtc` and  `output_traj` are actually used) and the output is made by concateneating the chains in ABCD order. Similarily for NMDA they are considered to be the same as the `input_pdb` are the same string, but they are read in as BCDA order.

To do the alignment separately different chains, we sould need separate input files for the laignment and fitting. For example in the above case, the AMPAR has four identical chains, while NMDAR has A and C chains matching, as well as B and D. As such we would normally want to do alignment of the two different NMDAR chains separately as we are only interested in their similarities to AMPAR, not between themselves. For the alignment we would have two input files, the first one as follows:
```
 input_pdb            | input_xtc            | selection | output_ndx              | output_traj

 input_files/ampa.pdb | input_files/ampa.xtc | chainID A | output_files/ampa_A.ndx | output_files/ampa_fitted.xtc
 input_files/ampa.pdb | input_files/ampa.xtc | chainID C | output_files/ampa_C.ndx | output_files/ampa_fitted.xtc

 input_files/nmda.pdb | input_files/nmda.xtc | chainID B | output_files/nmda_B.ndx | output_files/nmda_fitted.xtc
 input_files/nmda.pdb | input_files/nmda.xtc | chainID D | output_files/nmda_D.ndx | output_files/nmda_fitted.xtc
```

And similarily for the two other pairs of chains. For the fitting we could use the same input as before.