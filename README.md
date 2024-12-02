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
