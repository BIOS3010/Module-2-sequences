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
**X. .reverse() and .reverse_complement() ...**
- Compare with what we learned in class

**X. And so on ...**

**X. Working with FASTA format**

**6. Downloading some actual data (FASTA) and reading into Seq object**
