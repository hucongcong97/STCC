# STCC
consensus clustering enhances spatial domain detection for spatial transcriptomics data.

![image-20240213115617896](STCC/STCC.png)

# 1. Installation

```
pip install STCC==1.0.2
```

# 2. Data preparation

The input of STCC is a matrix composed of label vectors of different clustering results, where **the rows represent spots** and **the columns represent the indices of different clustering results**.

```
df = pd.read_csv('STCC/test data.csv',index_col=0)
```

![image-20240213124443246](STCC/data_display.png)

# 3. Running STCC

```
n_clusters = 15   # cluster numbers
method = 'wNMF-based' #  'Average-based','Onehot-based','wNMF-based'
seed = 2024 # random seed
if method == 'wNMF-based':
    labels_consensus, contributions = consensus_STCC(df,n_clusters,methods=method,seed=seed)
    print(f'Consensus labels are: {labels_consensus}')
    print(f'weight of {method} is:{contributions}')
else:
    labels_consensus = consensus_STCC(df,n_clusters,methods=method,seed=seed)
    print(f'Consensus labels are: {labels_consensus}')
```

![image-20240826204932456](STCC/results.png)

