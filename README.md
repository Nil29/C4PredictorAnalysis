<h1 align="left">C4PredictorAnalysis</h1>
Machine learning framework for the identification of sequence determinants of the convergent trait, C4 photosynthesis, in grasses
<br/>
<br/>
Provided here is the code accompanying the manuscript:
<br/>
<br/>
Yogadasan, N., Doxey, A.C., Chuong, S.D.X. (2023) A machine learning framework identifies plastid-encoded proteins harboring C3 and C4 distinguishing sequence    information. Genome Biology and Evolution (GBE)  - In Review
<br/>
<br/>
<h1 align="left">Requirements</h1>

Python v 3.8 <br/>
Biopython v 1.79 <br/>
Scikit-Learn v 1.1.1<br/>
Other dependencies as listed in .py files<br/>

<br/>

<h1 align="left">General Workflow</h1>
To reproduce summary statistics of v1, v2, and v3 models as reported in the manuscript...<br/>
<br/>
## v1 models (v1Modelling.py):<br/>
<br/>
Set the file_name variable to the correct alignment file you wish to query.
All alignment data is provided.<br/><br/>

```python
#example
file_name = 'alignment_data/rbcL_aligned.fasta' 
```
<br/>

Run v1Modelling.py <br/><br/>
Repeat for alignment data of interest

<br/>
v2 models (GetTenStrongestFeatures.py -> v2Modelling.py):<br/>
<br/>
First determine the 10 strongest features from v1 models. In GetTenStrongestFeatures.py
set the file_name variable to the correct alignment file you wish to query.

```python
#example
file_name = 'alignment_data/rbcL_aligned.fasta' 
```
<br/>

Run GetTenStrongestFeatures.py <br/><br/>

In v2Modelling.py set the file_name variable and the strongest features returned from GetTenStrongestFeatures.py

```python
#example
file_name = 'alignment_data/rbcL_aligned.fasta' 
strongest_features = [281, 468, 143, 328, 262, 228, 418, 101, 309, 270]
```

