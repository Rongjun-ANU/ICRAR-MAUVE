# 20260513 `ppxf_gas_wrapper.py` index error

The error occurs because the pipeline script `ppxf_gas_wrapper.py` has hardcoded variables that strictly assume there are exactly **three** "free" kinematic gas components. 

When our config sets `38 SIII` to `k17` (tied to line 17), we have 3 free components. Setting it to `f` (free) introduces a **4th** kinematic gas component.

Because of this mismatch, `pPXF` silently fails inside its worker routine on every spatial bin. It returns `NaN`s, and the pipeline later tries to fetch the 4th component's results (`component_x = 3` using zero-based indexing) from an array that was only built to store 3 size combinations, resulting in your `IndexError`.

That means if we want to use more than 3 free kinematic groups, we need to edit `ppxf_gas_wrapper.py` to make it dynamically size the arrays based on the number of gas components.

Around **Line 895**, it is where the hardcoding occurs:
```python
    n_gas_comp = 3  # len(np.unique(tpl_comp[gas_comp]))
```
And earlier around **Lines 849-861** and **875-888**, the starting kinematic guesses (`start` and `fixed` arrays) are visibly padded with exactly three hardcoded sets of gas parameters:
```python
            s = [
                ppxf_data[i][: config["KIN"]["MOM"]],
                [ppxf_data[i][0], 50],
                [ppxf_data[i][0], 50],
                [ppxf_data[i][0], 50],  # exactly 3 gas components
            ]
```

**The solution:**
Change the manual hardcoding to dynamically build the lists using `n_gas_comp`:

1. In `ppxf_gas_wrapper.py`, define the dynamic length before the loops by uncommenting the author's note (above line 845):
   ```python
   n_gas_comp = len(np.unique(tpl_comp[gas_comp]))
   ```
2. For the `start` list logic inside the `for i in range...` loops, replace it with a dynamic generator:
   ```python
   s = [ppxf_data[i][: config["KIN"]["MOM"]]] + [[ppxf_data[i][0], 50]] * n_gas_comp
   ```
3. For the `fixed` arrays:
   ```python
   f = [[1, 1, 1, 1]] + [[0, 0]] * n_gas_comp
   ```
4. Update **Line 895**:
   ```python
   # Remove the `n_gas_comp = 3` so it uses the dynamic variable you just declared.
   ```

