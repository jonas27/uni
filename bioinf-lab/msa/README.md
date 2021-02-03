# msa: Multiple Sequence Alignment

msa computes the alignments of three sequences. 

## Usage
*msa* can be only used as a function.:
```python
a1,a2,a3 = msa.run(seq1,seq2,seq3,cost_open,cost_extend,substitions)
```  

**Importantly**, all three sequences need to be passed as a char list with leading empty item.  
Example: `'AAT' should be ['','A','A','T']`  

## Tests
There is one test file `test_msa.py`. The test outputs the result for the assignment.  

To run the tests simple run a file:  
```bash
python3 test_msa.py
```