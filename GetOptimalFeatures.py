# Dependencies
from Bio import SeqIO

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


import numpy as np
from statistics import mode
import random
import sys
sys.setrecursionlimit(2000)


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


# Custom recursive feature elimination function (see manuscript methods)
def recursive_feature_elimination(X, y, length_remaining, strongest_indices, lvf=0, max_acc=0.5, remaining_residues = [], data={}):

    if length_remaining <= 2:
        print(strongest_indices, max_acc)
        return strongest_indices, max_acc, data
    
    new_strongest_indices = strongest_indices.copy()

    if remaining_residues == []:
        new_remaining_residues = [i+1 for i in range(length_remaining)]
        new_data = X
    else:
        del remaining_residues[lvf]
        new_remaining_residues = remaining_residues
        new_data = np.delete(X, lvf, 1)

    minimum_indices = []
    log_regression_scores = []
    

    log_classifier = LogisticRegression(max_iter = 300)
    scaler = StandardScaler()

    for k in range(501,600):
        x_train, x_test, y_train, y_test = train_test_split(new_data, y, test_size = 0.3, random_state = k)
        x_train_scaled = scaler.fit_transform(x_train)
        x_test_scaled = scaler.transform(x_test)

        log_classifier.fit(x_train_scaled, y_train)
        log_regression_scores.append(log_classifier.score(x_test_scaled, y_test))


        vector_coefficients = []
        for score in log_classifier.coef_.reshape(-1,1):
            vector_coefficients.append(abs(score[0]))
        minimum_index = vector_coefficients.index(min(vector_coefficients)) 
        minimum_indices.append(minimum_index)

    average_score = round((sum(log_regression_scores)/len(log_regression_scores)), 4)
    least_valuable_feature = mode(minimum_indices)
    data[len(new_remaining_residues)] = [average_score, new_remaining_residues.copy()]

    if average_score >= max_acc:
        max_acc = average_score
        new_strongest_indices = new_remaining_residues.copy()


    new_length_remaining = len(new_remaining_residues)

    print(len(new_strongest_indices), len(new_data[0]), max_acc, average_score)

    return recursive_feature_elimination(new_data, y, new_length_remaining, new_strongest_indices, least_valuable_feature, max_acc, new_remaining_residues, data)
    
    
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

c4_c3 = c4 + c3


c4_c3_vectors = c4_vectors + c3_vectors
data = np.array(c4_c3_vectors)


print('# C4 species:' + str(len(c4)))
print('# C3 species:' + str(len(c3)))
print('# of Initial Features:' + len(c4_vectors[0]))


# for randomization
rand_c4_c3 = c4_c3.copy()
def randFunction():
    return 0.5
random.shuffle(rand_c4_c3, randFunction)


strongest, score, total_data = recursive_feature_elimination(data, c4_c3, len(c4_vectors[0]), [])
print(strongest, score, total_data[10])

# for permutation analysis comparison - use rand_c4_c3 as argument instead of c4_c3









