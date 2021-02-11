# Exercise 2: Using Biopython

In this exercise, you will use what you learned in the [introductory exercise](Exercise1.md) to analyse some real data.

## 2.4.1 Downloading some actual data (FASTA) and read it into a Seq object
- Go to NCBI - Gene: https://www.ncbi.nlm.nih.gov/gene/
- Search for the human beta globin gene (HBB)
- Download the HBB gene in FASTA format and save it on your computer as "sequence.fasta"
- Load this FASTA file into a Seq object called `hbb`. This can be done by using the SeqIO module in Biopython. This module has a function `read` that takes two arguments, a file name and a file type (in our case "fasta"):
```python
from Bio import SeqIO
seq_record = SeqIO.read("sequence.fasta", "fasta")
# The actual sequence is in the seq variable of the seq_record object:
hbb = seq_record.seq
```
```diff
! Translate the sequence into the corresponding amino acid sequence. 
! Compare with the amino acid sequence on NCBI, is this as expected?
! Advanced: Use python string slicing on `hbb` to explore the translation further
```
## 2.4.2 Reading a FASTA-file containing multiple sequences into Seq objects

Start by downloading the coding DNA sequence of all human genes in FASTA format:

- Go to the UCSC Table Browser: https://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=250654185_KikUKaMDmiGmkTgDVb2wBpVbjc9K&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=knownGene&hgta_table=knownGene&hgta_regionType=genome&position=chrX%3A15%2C578%2C261-15%2C621%2C068&hgta_outputType=sequence&hgta_outFileName=all_genes.fasta
- Inspect the options that are selected
- Click `get output`
- Select `Genomic` and click Submit`
- Select only `CDS Exons` ("5' UTR Exons", "3' UTR Exons", and "Introns" should be unselected)
- Click `get sequence"`
- Save the `all_genes.fasta` file in the current working directory (the `Module-2-sequences directory`) (! Warning: the file is 105.3 MB)

To loop through all the genes and print the sequence (Seq object) do this:

```python
from Bio import SeqIO

all_genes = SeqIO.parse("all_genes.fasta", "fasta")

for gene in all_genes:
    print(gene.id)
    print(gene.seq)
```

```diff
+ Note:
+ When importing libraries in Python, we usually put the import statement at the top of the file. 
+ This makes what we have imported available in the whole file
+ The statement `from Bio.Seq import Seq` makes the `Seq` class available to us. 
+ An alternative way of importing is `import Bio`. This makes the `Seq`-class available through `Bio.Seq`.
```

## 2.4.3 Using a text editor to create Python scripts
To work efficiently with Python, it is good to create scripts that store our commands and allow us to run them in the Terminal. A python script is simply a text file containing the commands/functions we need to run in the right order. When our python code becomes multiple lines, for-loops, and so on, it is better to store these as files. This will will do now.

First, exit the interactive python session by using:
```bash
exit()
```

You may already be used to a particular text editor to create scripts. Feel free to use this, but if you are new to this you should use the [nano text editor](https://www.howtogeek.com/howto/42980/the-beginners-guide-to-nano-the-linux-command-line-text-editor/). This is a very simple text editor running inside our terminal.
To open it, type the following into your Terminal (in Bash):
```bash
nano
```
- Paste in the python code from the previous step (**2.4.2**) into the text editor
- Save the file (CTRL - Shift - O), choose filename 'printgenes.py'
- Exit `nano` (CTRL - Shift - X)

```diff
! Use what we learned the last week to verify that `printgenes.py` contains the python script
```
## 2.4.4 Running a python script from the terminal

Run the python script like this:
```python
python printgenes.py
```

```diff
! Use what we learned the last week to redirect the output into a file called `genes.txt`
```

## 2.4.5 Editing and expanding the python script
Modify the Python script by again opening it in your text editor (e.g. nano) and then:
- Expand/change the code in `printgenes.py`, so that it outputs the gene-ID and the protein length on the same line. You can skip printing the sequence in order to make the output more readable

```diff
! Use what we learned the last week to find the id of the longest gene
! And the shortest gene
```

## 2.4.6 Make a small tool that locates genes containing a specific sequence**
We now want to make a small Python tool that lets the user search for genes that contain a spescific DNA sequence (pattern). For convenience, we want to be able to run this tool on the command line like this:

```bash
python3 find_genes.py ACTG
```

- Make a new file `find_genes.py`  where you import SeqIo, read in all the genes and loop through them in the same way as in the previous exercise.
- Add the following line at the top of the file `import sys`. This gives us access to `sys.argv` which is a list containing the arguments specified on the command line when running the python script.
- You can read the first argument like this:
```python
pattern = sys.argv[1]
```

```diff
! Modify your code so that only the id of genes containing the pattern specified by the user are printed
! What happens if you run the program without specifying a pattern? 
! Can you modify your code so that it doesn't crash but instead prints an error message to the user?
! If no genes contain the pattern, print a descriptive message (e.g. "No genes found")
```
