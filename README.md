# Module-2-sequences
**1. Install Biopython**
Show here

**2. Load Biopython**

```python
import Bio
from Bio.Seq import Seq
```

**3. Create Seq object**
```python
my_seq = Seq("AGTACACTGGT")
```
```diff
! Test the python function `len` on the `my_seq`, is it working in the same way as before?
! What is the purpose of a Seq object? Why not just use a string?
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



