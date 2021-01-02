# Module-2-sequences
In this exercise, you will get to know how to represent and operate on DNA sequences using Biopython. In the [next exercise](Exercise2.md)
you will use what you learn here on a bigger problem.

#### Representing a DNA sequence as a Python string vs. using Biopython
You've previously learned that representing a DNA sequence in Python can be done using strings:
```python
my_seq = "ACTG"
```
... and that Python makes it easy to perform common operations on such strings, such as counting number of A's inside the DNA sequence, or getting the length of the sequence:
```python
length_of_sequence = len(my_seq)
number_of_a = my_seq.count("A")
```

However, if you want to do more complicated things, like getting the reverse complement of the sequence, or translating it to amino acids, that is not straight forward and would require some code. Luckily, the Biopython package supports all this functionality and much more, making  


**1. Import Biopython**

Make a file called `biopython_test.py` where you will try out Biopython.

In the top of the file, start by importing the `Seq` class from Biopython:

```python
from Bio.Seq import Seq
```

A few notes:
* When importing stuff in Python, we usually put the import statement at the top of the file. This makes what we have imported available in the whole file, and it makes it easy to get an overview of everything that is imported when all import statements are togheter at top of the file.
* The statement `from Bio.seq import Seq` makes the `Seq` class available to us. An alternative way of importing, that you may have seen is `import Bio`. That makes the `Seq`-class available through `Bio.Seq`.

**3. Create Seq object**

Use the following code to create a new sequence object:
```python
my_seq = Seq("AGTACACTGGT")
```
```diff
! Test the python function `len` on the `my_seq`, is it working in the same way as before?
! What do you think is the purpose of a Seq object? Why not just use a string?
! Try to see if you manage to convert the Seq-object to a string (without reading any documentation). 
```


**X. .complement() and .reverse_complement() ...**
```python
my_seq.complement()
my_seq.reverse_complement()
my_seq2 = Seq("AGYCCTD")
my_seq2.complement()
```
```diff
! Explain what what happens with `my_seq` and `my_seq2`
! Refer back to the lecture slides if you have forgotten about (reverse) complementarity
```

- Compare with what we learned in class

**X. .find() to search for sequences**
TBD.
```python
my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
my_rna.find("AUG")
my_rna.find("UAG")
my_rna.find("UAA")
my_rna.find("UGA")
```

**X. .count() to count occurences of sequences**
TBD.
```python
my_rna.count("AUG")
```
```diff
! Make a Seq object with the sequence "AAAA"
! Count the occurence of "AA" in this object
! Explain what happened
```

**X. Test the .transcribe() method**
```python
my_seq.transcribe()
```
```diff
! What is the difference between a method and a function?
! Explain what the `.transcribe()` method above does on a Seq object
! What happens if you try to transcribe `my_rna`?
```

**X. Test the .translate() method**
```python
my_rna.translate()
```
```diff
! Explain what the `.translate()` method above does on a Seq object
! Explain what happens if you translate the `my_seq`
```

**X. Working with FASTA format with SeqIO**

**6. Downloading some actual data (FASTA) and reading into Seq object**
- Go to NCBI - Gene: https://www.ncbi.nlm.nih.gov/gene/
- Search for the human beta globin gene (HBB)
- Download the HBB gene in FASTA format and save it on your computer as "sequence.fasta"
- Load this FASTA file into a Seq object called `hbb`:
```python
from Bio import SeqIO
seq_record = SeqIO.read("sequence.fasta", "fasta")
hbb = seq_record.seq
```
```diff
! Translate the sequence into the corresponding amino acid sequence. 
! Compare with the amino acid sequence on NCBI, is this as expected?
! Advanced: Use python string slicing on `hbb` to explore the translation further
```
**7. Reading a FASTA-file containing multiple sequences into Seq objects**
Start by download the coding DNA sequence of all human genes in FASTA format.

- Go to the UCSC Table Browser: https://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=250654185_KikUKaMDmiGmkTgDVb2wBpVbjc9K&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=knownGene&hgta_table=knownGene&hgta_regionType=genome&position=chrX%3A15%2C578%2C261-15%2C621%2C068&hgta_outputType=sequence&hgta_outFileName=all_genes.fasta
- Inspect the options that are selected
- Click `get output`
- Select `Genomic` and click Submit`
- Select only `CDS Exons` ("5' UTR Exons", "3' UTR Exons", and "Introns" should be unselected)
- Click `get sequence"`
- Save the `all_genes.fasta` file in your current working directory (Warning: the file is 105.3 MB)

To loop through all the genes and print the sequence (Seq object), one can do this:

```python
all_genes = SeqIO.parse("genes_all.fa", "fasta")

for gene in all_genes:
    print(gene.seq)    
```

```diff
! Print the protein length (number of amino acids) of each gene
! Calculate the average protein length (number of amino acids)
```



