### Install Biopython

The easiest way to install Biopython is to use the builtin Python package manageer (PiP).

If you have Python installed on your system, you usually also have PiP installed as a module, and can install Biopython executing this command in the terminal:

## 2.1.1 Opening the terminal
In the exercises last week, Windows-users installed the Gitbash terminal. Mac (and Linux) users already have a Terminal application installed in their systems.

- Open your Terminal application (Gitbash for Windows users)

## 2.1.2 Installing Biopython using pip
Execute the following command in the terminal:
```
python3 -m pip install biopython
```

```diff
+ Note: 
+ You may need to swap `python3` with `python` in the command above, if you get an error.
```

## 2.1.3 Opening an interacive python session
We can interact with Python in our Terminal, just as we do with Bash (and in Python notebook).
Type `python3` (or `python`) on the terminal to go into an interactive python shell:

```bash
python3
```

```diff
+ Note:
+ the >>> symbols to the left, indicate that we are interacting with Python
```

## 2.1.4 Testing the Biopython installation
Type the following into the terminal
```python
import Bio
```
If no errors are shown, Bipython is installed. 

## 2.1.5 Interacting with Python in the terminal
Try out some simple Python commands, once you are inside the interactive Python session:
```python
1+1
a = 1
print("hei")
```
Feel free to try other commands as well

```diff
! What is the difference between interacting with python and running commands in Bash?
```

## 2.1.6 Exiting the interactive Python session:
Type `exit()` to leave the interactive python shell:
```python
exit()
```
```diff
+ Note that you are now back in the shell, and can run Bash commands again
```



