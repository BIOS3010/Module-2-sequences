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

**4. Test the .transcribe() method**
```python
my_seq.transcribe()
```
```diff
! Explain what the `.transcribe()` method above does on a Seq object
! What is the difference between a method and a function?
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
```
**X. .count() to count occurences of sequences**
TBD.
```python
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
TBD.
```python
```
```diff
! Loop through sequences and search for start codons
! Loop through sequences and count the occurence of 'GC'
```



