# Gotoh

Gotoh is a python library for computing the score and alignments of two sequences. 

## Usage
*Gotoh* can be used as a function or class.  

As a class:
```python
import gotoh

g = gotoh.Gotoh()
score, alignments = g.run(sequence_1,sequence_2,t[2],t[3],cap_cost_open, cap_cost_extend, substition_matrix)
```
As a function:
```
g = gotoh.gotoh()
```  

**Importantly**, both sequences need to be passed as a char list with leading empty item.  
Example: `'AAT' should be ['','A','A','T']`  
However in the file *gotoh_helpers.py* you can find a function for parsing fasta files into the expected format.  
Example:
```python
from gotoh_helpers import read_fasta_file

sequence = read_fasta_file('<file_path>')
```

## Assumptions
The algorithm expects a minus inifinty min score for the initialization of P and Q.
In order to avoid imports, a theoretical min score of -1,000,000, -10^6 was set.
Also, in order to avoid further checks for `None` checks for the same matrices
it uses a score of -10^10 to mark never reached fields.

## Tests
There are two test fiels `test_gotoh.py` and `test_gotoh_unit.py`, respectively. The former tests the full gotoh algorithm while the latter tests single functions from both, the `gotoh.py` and `gotoh_helpers.py`.  
To run the tests simple run a file:  
```bash
python test_gotoh.py && \
python test_gotoh_unit.py 
```