# Dependencies
from Bio import SeqIO

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

import numpy as np
import random
from statistics import pstdev

# Filename with the relevant alignment data
# Ensure that sequence descriptions in fasta file have appropriate labels included (ex. ">C4_MPC_Paspalidium_geminatum" "C4_" or "C3_" in description will be used to organize data later)
file_name = 'rbcL_aligned.fasta' 

# Functions for feature representation of alignment data
def create_sequence_subsets(sequence):
    sequence_split = [sequence[x:x+1] for x in range(0, len(sequence), 1)]
    return sequence_split

    # Feature engineering function (see manuscript methods)
def create_map(split_sequences):
    position_maps = []
    for i in range(len(split_sequences[0])):
        map = {}
        value = 0
        for sequence in split_sequences:
            if sequence[i] not in map:
                map[sequence[i]] = value
                value += 1
        position_maps.append(map)
    return position_maps 

def create_vector(split_sequence, position_map):
    vector = []
    for i in range(len(split_sequence)):
        map = position_map[i]
        vector.append(map[split_sequence[i]])
    return vector 


# Converting sequences in alignment to feature vector representations
sequences = []
descriptions = []
for seq_record in SeqIO.parse(file_name, 'fasta'):
    sequences.append(seq_record.seq)
    descriptions.append(seq_record.description)

split_sequences = []
for sequence in sequences:
    split_sequences.append(create_sequence_subsets(sequence))

feature_map = create_map(split_sequences)

split_sequences_vectors = []
for sequence in split_sequences:
    split_sequences_vectors.append(create_vector(sequence, feature_map))


# Organize sequences from alignment according to classification label in descriptions
c4_vectors = []
c3_vectors = []

    # If "C4_" in sequence description, append vector to C4 set of vectors
    # If "C3_" in sequence description, append vector to C3 set of vectors
for j in range(len(split_sequences_vectors)):
    if 'C4_' in descriptions[j]:
        c4_vectors.append(split_sequences_vectors[j])
    if 'C3_' in descriptions[j]:
        c3_vectors.append(split_sequences_vectors[j])


# Binary classification labels - C3: 0, C4: 1
c4 = [1 for i in range(len(c4_vectors))]
c3 = [0 for i in range(len(c3_vectors))] 


print('# C4 species:' + str(len(c4)))
print('# C3 species:' + str(len(c3)))
print('# of Initial Features:' + len(c4_vectors[0]))

# Y-data for train/test
c4_c3 = c4 + c3

# x-data for train/test
c4_c3_vectors = c4_vectors + c3_vectors
data = np.array(c4_c3_vectors)


# Initial modelling to return 10 strongest features
log_classifier = LogisticRegression(max_iter = 300)
scaler = StandardScaler()


# Minimal Test Set for initial modelling so that 10 strongest features can be inferred from all sequence data
x_train, x_test, y_train, y_test = train_test_split(data, c4_c3, test_size = 0.001)
x_train_scaled = scaler.fit_transform(x_train)
x_test_scaled = scaler.transform(x_test)

log_classifier.fit(x_train_scaled, y_train)


# Organizing regression coefficients by absolute value to determine 10 strongest features (greatest coefficient values)
vector_coefficients = []

for score in log_classifier.coef_.reshape(-1,1):
    vector_coefficients.append(abs(score[0]))

coefficients = []
i = 1
for score in vector_coefficients:
    coefficients.append([score, i])
    i+= 1

coefficients.sort(reverse = True)

ymax1, xmax1 = coefficients[0][0], coefficients[0][1]
ymax2, xmax2 = coefficients[1][0], coefficients[1][1]
ymax3, xmax3 = coefficients[2][0], coefficients[2][1]
ymax4, xmax4 = coefficients[3][0], coefficients[3][1]
ymax5, xmax5 = coefficients[4][0], coefficients[4][1]
ymax6, xmax6 = coefficients[5][0], coefficients[5][1]
ymax7, xmax7 = coefficients[6][0], coefficients[6][1]
ymax8, xmax8 = coefficients[7][0], coefficients[7][1]
ymax9, xmax9 = coefficients[8][0], coefficients[8][1]
ymax10, xmax10 = coefficients[9][0], coefficients[9][1]

ymax = [ymax1, ymax2, ymax3, ymax4, ymax5, ymax6, ymax7, ymax8, ymax9, ymax10]

xmax = [xmax1, xmax2, xmax3, xmax4, xmax5, xmax6, xmax7, xmax8, xmax9, xmax10]

print('10 Strongest Features: ' + xmax)




