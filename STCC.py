import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from scipy import stats
from cvxopt import matrix, solvers
import anndata
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances as pair
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import AgglomerativeClustering
from sklearn.manifold import SpectralEmbedding
from scipy.stats import pearsonr,spearmanr
from sklearn.decomposition import PCA,NMF
from numpy.linalg import  *

import math
import time
import random
import os

df = pd.read_csv('test data.csv',index_col=0)
n_clusters = 12   # true number of categories
method = 'Averaged-based' # 'Onehot-based','NMF-based','wNMF-based'
seed = 2024 # random seed
labels_consensus = consensus_summary(df,n_clusters,methods=method,seed=seed)