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

# Input strongest features for alignment data (as a list) returned from GetTenStrongestFeatures.py 
strongest_features = []

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

    # Function to convert full-length vectors to vectors of only strongest features, given strongest feature positions
def create_strongest_feature_vector(split_vector, features):
    strongest_vector = []
    for position in features:
        strongest_vector.append(split_vector[position - 1])
    return strongest_vector


# Converting sequences in alignment to strongest feature vector representations
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

strongest_feature_vectors = []
for sequence in split_sequences_vectors:
    strongest_feature_vectors.append(create_strongest_feature_vector(sequence, strongest_features))

# Organize sequences from alignment according to classification label in descriptions
c4_vectors = []
c3_vectors = []

    # If "C4_" in sequence description, append vector to C4 set of vectors
    # If "C3_" in sequence description, append vector to C3 set of vectors
for j in range(len(split_sequences_vectors)):
    if 'C4_' in descriptions[j]:
        c4_vectors.append(strongest_feature_vectors[j])
    if 'C3_' in descriptions[j]:
        c3_vectors.append(strongest_feature_vectors[j])


# Binary classification labels - C3: 0, C4: 1
c4 = [1 for i in range(len(c4_vectors))]
c3 = [0 for i in range(len(c3_vectors))] 

c4_c3 = c4 + c3


c4_c3_vectors = c4_vectors + c3_vectors
data = np.array(c4_c3_vectors)

print('# C4 species:' + str(len(c4)))
print('# C3 species:' + str(len(c3)))
print('# of Features:' + len(c4_vectors[0]))



# For logistic regression model building and repeated random sub-sampling (see manuscript methods)

log_classifier = LogisticRegression(max_iter = 300)
scaler = StandardScaler()


log_regression_scores = []
fprs = []
tprs = []
auc_scores = []

# Generating test scores for 500 different and random splits of data into Train/Test, using test size of 30%
for k in range(501,1000):
    x_train, x_test, y_train, y_test = train_test_split(data, c4_c3, test_size = 0.3, random_state = k)
    x_train_scaled = scaler.fit_transform(x_train)
    x_test_scaled = scaler.transform(x_test)
    log_classifier.fit(x_train_scaled, y_train)
    log_regression_scores.append(log_classifier.score(x_test_scaled, y_test))

    y_pred = log_classifier.predict_proba(x_test_scaled)
    y_pred = y_pred[:,1]
    fpr, tpr, _ = roc_curve(y_test, y_pred)
    auc = roc_auc_score(y_test, y_pred)
    fprs.append(fpr)
    tprs.append(tpr)
    auc_scores.append(auc)\
    
average_score = round((sum(log_regression_scores)/len(log_regression_scores)), 4)
stdev_acc = round(pstdev(log_regression_scores),4)
average_auc = round((sum(auc_scores)/len(auc_scores)),4)
stdev_auc = round(pstdev(auc_scores),4)


# For random permutation comparison
rand_c4_c3 = c4_c3.copy()
random.shuffle(rand_c4_c3)

log_regression_scores_rand = []
fprs_rand = []
tprs_rand = []
auc_scores_rand = []

for h in range(501,1000):
    x_train_rand, x_test_rand, y_train_rand, y_test_rand = train_test_split(data, rand_c4_c3, test_size = 0.3, random_state = h)
    x_train_rand_scaled = scaler.fit_transform(x_train_rand)
    x_test_rand_scaled = scaler.transform(x_test_rand)
    log_classifier.fit(x_train_rand_scaled, y_train_rand)
    log_regression_scores_rand.append(log_classifier.score(x_test_scaled, y_test))

    y_pred_rand = log_classifier.predict_proba(x_test_rand_scaled)
    y_pred_rand = y_pred_rand[:, 1]
    fpr_rand, tpr_rand, _ = roc_curve(y_test_rand, y_pred_rand)
    auc_rand = roc_auc_score(y_test_rand, y_pred_rand)
    fprs_rand.append(fpr_rand)
    tprs_rand.append(tpr_rand)
    auc_scores_rand.append(auc_rand)

average_score_rand = round((sum(log_regression_scores_rand)/len(log_regression_scores_rand)), 4)
stdev_acc_rand = round(pstdev(log_regression_scores_rand),4)
average_auc_rand = round((sum(auc_scores_rand)/len(auc_scores_rand)), 4)
stdev_auc_rand = round(pstdev(auc_scores_rand),4)


print('Average Classification Accuracy: ' + average_score)
print('Average Classification Accuracy (permutation): ' + average_score_rand)




