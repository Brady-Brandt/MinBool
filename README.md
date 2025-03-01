## Example 
```python3
from minbool import get_expr, get_truth_table
minterms = [0,1,2,3]
dont_care = [7]
# number of input variables
num_vars = 3

get_truth_table(num_vars, minterms, dont_care)
print("F = ", get_expr(num_vars, minterms, dont_care))
```

```sh
Dec | A B C  O/P
-----------------
0   | 0 0 0 | 1
1   | 0 0 1 | 1
2   | 0 1 0 | 1
3   | 0 1 1 | 1
4   | 1 0 0 | 0
5   | 1 0 1 | 0
6   | 1 1 0 | 0
7   | 1 1 1 | X
F =  A'
```

## Support for Finite State Machines

### Add your own State Table 
- Dictionary where each key represents the current state 
- Value is a list of strings representing the next state and output  
- Index 0 & 1 denote the next state when x=0 or x=1  
- For a mealy machine index 2 & 3 denote the output when x = 0 or x = 1 
- Don't cares are represented by an x
- Dictionary keys can only be integers the backend will take care of the don't cares

### Moore Machine
  
  ```python3
from minbool import * 
fsm = FiniteStateMachine(is_mealy=False)
state = {
    0: ['0','1','0'],
    1: ['0','2','0'],
    2: ['3','2','0'],
    3: ['0','4','0'],
    4: ['0','2','1'],
}

fsm.add_state_table(state)
fsm.get_boolean_expressions()
```
#### Output
```sh
A+ = BCX
B+ = BC' + AX + B'CX
C+ = BC'X' + A'B'C'X
Z =  A
```

### Mealy Machine 
```python3
from minbool import * 

fsm = FiniteStateMachine(is_mealy=True)

state = {
    0: ['0','1','0', '0'],
    1: ['2','1','0', '0'],
    2: ['0','1','0', '1'],
}

fsm.add_state_table(state)
fsm.get_boolean_expressions()
```
#### Output
```sh
A+ = BX'
B+ = X
Z =  AX
```

## Sequence Detectors 
Sequence detectors are a subclass of FiniteStateMachine 

```python3
from minbool import *

seq = SequenceDetector('1101', is_mealy=False)
seq.print_state_table()
seq.get_boolean_expressions()
```

#### Output
```sh
S0    S0    S1    0
S1    S0    S2    0
S2    S3    S2    0
S3    S0    S4    0
S4    S0    S2    1
A+ = BCX
B+ = BC' + AX + B'CX
C+ = BC'X' + A'B'C'X
Z =  A
```





