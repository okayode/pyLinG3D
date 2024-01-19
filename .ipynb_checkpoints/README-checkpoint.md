# pyLinG3D: Visualizing the Spatio-Temporal Dynamics of Clonal Evolution

This library is a python implementation of the LinG3D (https://github.com/rejniaklab/LinG3D) which generates the 3D lineage trees of (1) all clones; (2) individual clones; (3) all clones, but with only those cells that survived to the end of simulation; and (4) individual clones containing only those cells that survived to the end of simulation.

## Installing with pip

if you are using Linux or macOS you can install pyLinG3D with pip

```bash
pip install pyLinG3D
```
## Example

```python
%matplotlib ipympl

from pyLinG3D import LinG3DClone
from platform import python_version
print(python_version())
```

```python
#exampleB05

pathData='exampleB05'
cloneNum=5
toPrint=1                # save the final figure
IsGradient = 1           # draw drug gradient in the background 1-yes; 0-no;
xmin=-100; xmax=100; ymin=xmin; ymax=xmax  # 2D domain boundaries
tmin=0; tmax=100000                        # time/iteration boundaries
fileStep = 2000       # frequency of data 
    
LinG3DClone(pathData,cloneNum,IsGradient,xmin,xmax,ymin,ymax,tmin,tmax,fileStep,toPrint).linG3DClone()
```
<div style="margin:2%";>  
    <img src="exampleB05\fig_clones\tree_clone_5.jpg"; alt="tree_clone_5"; width=50%;/>
</div>

### Authors
Anjun Hu, Maureiq Ojwang’, Kayode Olumoyin kayode.olumoyin@moffitt.org, Katarzyna Rejniak 

### Source Code
https://github.com/okayode/pyLinG3D


### Acknowledgements

This work was supported by the NIH/NCI grants CA259387 and CA272601, and the NIH/NCI Physical Sciences Oncology Network (PSON) grant CA202229. 

