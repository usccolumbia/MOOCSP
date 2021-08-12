# MOOCSP: age-fitness multiobjective GA for crystal structure prediction

Introduction\\
The code provides 3 major functions:

Generate the input file from the cif file

Crystal Structure prediction by contact map,the coordination number of the cations and the ages of the individuals


Evaluate the predicted structure,calculate contact map accuracy,coordination error,rmsd,mae and rms

Installation
Install the relevant packages if not already installed:

pymoo (0.4.2.2)
numpy (1.17.2)
scikit-learn (0.21.3)
pyxtal(0.2.2)
pytmatgen (2020.11.11)

pymoo, numpy, scikit-learn, pyxtal, and pymatgen
pip install pymoo
pip install numpy
pip install scikit-learn
pip install pyxtal
pip install pymatgen

Since we have improved the pymoo algorithm, we need to overwrite some source files in the pymoo library.
python cover.py

replace follow files：
D:\lib\site-packages\pymoo\algorithms\nsga2.py
D:\lib\site-packages\pymoo\algorithms\genetic_algorithm.py
D:\lib\site-packages\pymoo\model\mating.py
D:\lib\site-packages\pymoo\model\infill.py


Dataset
The data that support the findings of this study are openly available in Materials Project database at http:www.materialsproject.org


Usage

Generate the input file from the cif file
python cm_test.py --cif 2-14-mp-236.cif

Crystal Structure prediction by contact map,the coordination number of the cations and the ages of the individuals
python age_fitness.py --input 2-14-mp-236.input

Evaluate the predicted structure,calculate contact map accuracy,coordination error,rmsd,mae and rms
python measure.py --cif 2-14-mp-236.cif --predicted 2-14-mp-236_predicted.cif


or python exp1_batch.py
