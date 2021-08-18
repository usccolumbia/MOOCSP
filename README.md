## MOOCSP: age-fitness multiobjective GA for crystal structure prediction

Source code for our paper "Crystal structure prediction using age-fitness multi-objective genetic algorithm and coordination number constraints".
https://arxiv.org/abs/2107.01346

Wenhui Yang (Guizhou University)\
Jianjun Hu (University of South Carolina, Machine learning and Evolution Laboratory)

### Introduction

The code provides 3 major functions:

Generate the input file from the cif file

Crystal Structure prediction by contact map,the coordination number of the cations and the ages of the individuals. The Genetic algorithm can be run in parallel with multi-core CPU.


Evaluate the predicted structure,calculate contact map accuracy,coordination error,rmsd,mae and rms

### Installation

Step1: create a virtual environment
```console
python -m venv py3
source py3/bin/activate
```
Step2:\
Install the relevant packages if not already installed:
```console
pip install pymoo numpy scikit-learn pyxtal and pymatgen
```
pymoo (0.4.2.2)\
numpy (1.17.2)\
scikit-learn (0.21.3)\
pyxtal(0.2.2)\
pytmatgen (2020.11.11)

or individually:
```console
pip install pymoo
pip install numpy
pip install scikit-learn
pip install pyxtal
pip install pymatgen
```
Step3: \
Since we have improved the pymoo GA algorithm, we need to overwrite some source files in the pymoo library.\
First, run  
```console
pip -V  to figure out the python package path for pymoo:\
```
If your output is: pip 20.2.4 from /Users/xxx/opt/anaconda3/lib/python3.7/site-packages/pip (python 3.7)\
then, the pymoo path is /Users/xxx/opt/anaconda3/lib/python3.7/site-packages/pymoo

replace the following files by running the following commands:\

```console
export PYMOODIR=/Users/xxx/opt/anaconda3/lib/python3.7/site-packages/pymoo
cp nsga2.py ${PYMOODIR}/algorithms/nsga2.py
cp genetic_algorithm.py ${PYMOODIR}/algorithms\genetic_algorithm.py
cp mating.py ${PYMOODIR}/model/mating.py
cp infill.py ${PYMOODIR}/model/infill.py
```

<!-- python cover.py -->
or run:
```console
patch.sh /Users/xxx/opt/anaconda3/lib/python3.7/site-packages/pymoo
```
### Dataset
The data that support the findings of this study are openly available in Materials Project database at http:www.materialsproject.org


### Usage

Generate the input file from the cif file\
$python cm_test.py --cif 2-14-mp-236.cif

Crystal Structure prediction by contact map,the coordination number of the cations and the ages of the individuals\
$python age_fitness.py --input 2-14-mp-236.input

Evaluate the predicted structure,calculate contact map accuracy,coordination error,rmsd,mae and rms\
$python measure.py --cif 2-14-mp-236.cif --predicted 2-14-mp-236_predicted.cif


or python exp1_batch.py
