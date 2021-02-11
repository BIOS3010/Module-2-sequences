## Exercise 2: Using Biopython

In this exercise, you will use what you learned in the [introductory exercise](Exercise1.md) to analyse some real exercises. In the last exercise you will make a small tool that automatically fetches sequences from .


**1. Downloading some actual data (FASTA) and read it into a Seq object**
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
**2. Reading a FASTA-file containing multiple sequences into Seq objects**

Start by download the coding DNA sequence of all human genes in FASTA format:

- Go to the UCSC Table Browser: https://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=250654185_KikUKaMDmiGmkTgDVb2wBpVbjc9K&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=knownGene&hgta_table=knownGene&hgta_regionType=genome&position=chrX%3A15%2C578%2C261-15%2C621%2C068&hgta_outputType=sequence&hgta_outFileName=all_genes.fasta
- Inspect the options that are selected
- Click `get output`
- Select `Genomic` and click Submit`
- Select only `CDS Exons` ("5' UTR Exons", "3' UTR Exons", and "Introns" should be unselected)
- Click `get sequence"`
- Save the `all_genes.fasta` file in your current working directory (Warning: the file is 105.3 MB)

To loop through all the genes and print the sequence (Seq object), one can do this:

```python
from Bio import SeqIO

all_genes = SeqIO.parse("all_genes.fasta", "fasta")

for gene in all_genes:
    print(gene.id)
    print(gene.seq)
```

```diff
! Expand the code above so that also the protein length (number of amino acids) of each gene is printed (you can skip printing the sequence in order to make the output more readable)
! Calculate the average protein length (number of amino acids) of each gene
```

**3 Make a small tool that locates genes containing a specific sequence**
We now want to make a small Python tool that lets the user search for genes that contain a spescific DNA sequence (pattern). For convenience, we want to be able to run this tool on the command line like this:

```bash
python3 find_genes.py ACTG
```

* Make a new file `find_genes.py`  where you import SeqIo, read in all the genes and loop through them in the same way as in the previous exercise.
* Add the following line at the top of the file `import sys`. This gives us access to `sys.argv` which is a list containing the arguments specified on the command line when running the python script.
* You can read the first argument like this:
```python
pattern = sys.argv[1]
```

```diff
! Modify your code so that only the id of genes containing the pattern specified by the user are printed
! What happens if you run the program without specifying a pattern? Can you modify your code so that it doesn't crash but instead prints an error message to the user?
! Also make sure that the tool prints the average length of the genes that contain the pattern. If no genes contain the pattern, print a descriptive message (e.g. "No genes found")
```
