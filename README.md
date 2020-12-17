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

**X. And so on ...**

**X. Working with FASTA format**

**6. Downloading some actual data (FASTA) and reading into Seq object**
