# Working with DNA sequences using Biopython
In this exercise, you will get to know how to represent and operate on DNA sequences using Biopython. In the next exerise
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

However, if you want to do more complicated things, like getting the reverse complement of the sequence, or translating it to amino acids, that is not straight forward and would require some code. Luckily, the Biopython package supports all this functionality and much more.


## 2.3.1 Opening an interactive Python session
As we did in the previous exercise, open Python interactively by running this in your Terminal:
```python
python3
```
## 2.3.2 Importing the Seq class from Biopython
Type this in
```python
from Bio.Seq import Seq
```

```diff
+ Note:
+ The statement `from Bio.Seq import Seq` makes the `Seq` class available to us. 
+ An alternative way of importing is `import Bio`. This makes the `Seq`-class available through `Bio.Seq`.
```

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


**3. .complement() and .reverse_complement() ...**

By calling the `complement()` or `reverse_complement()` methods on the Seq-object, we can get the complement or reverse complement of a sequence:
```python
my_seq_complement = my_seq.complement()
my_seq_reverse_complement = my_seq.reverse_complement()
```
```diff
! Run the code above and also print the two new sequence `my_seq_complement` and `my_seq_reverse_complement`. Is the result as expected?
! Is the original variable `my_seq` changed in any way? 
! Refer back to the lecture slides if you have forgotten about (reverse) complementarity
```

- Compare with what we learned in class

**4. .find() to search for sequences**

Bipython also makes it easy to find a subsequence within a sequence by using the `find`-method. This method takes one argument (the subsequence to be found).

```diff
! Defining a new Biopython sequence (you can choose what the sequence should be)
! Use the `find()` method to find the location of a subsequence. Remember that you use methods on the sequence object by using the `some_sequence.some_method()`-notation. Store the result in a variable and print that variable to see the location.
! What happens if the subsequence you try to find does not exist in the sequence?
! What happens if the subsequence exists multiple times?
```


**5. Other methods**

In addition to the methods introduced above, the Seq-object has other methods that can be used:
* `sequence.count(subsequence)`: Returns the number of times subsequence occurs in the sequence
* `sequence.transcribe()`: Transcribes the sequence from DNA til RNA
* `sequence.translate()`: Translates from DNA to RNA

```diff
! Play around with these methods to see how they work
! What happens if you try to transcribe an RNA sequence or translate a DNA sequence?
```

**X. Working with FASTA format with SeqIO**


Go to the [next exercise](Exercise2.md)
