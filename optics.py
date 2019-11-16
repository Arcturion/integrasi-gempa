from sklearn.cluster import OPTICS, cluster_optics_dbscan
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wradlib as wrl

# Generate sample data


file_radar = wradlib.util.get_wradlib_data_file('')
img, metadata = wrl.io.read_GAMIC_hdf5(file_radar)
ax1, pm = wradlib.vis.plot_ppi(img)

data_balai = pd.read_csv('balai.csv', sep=',')
data_dramaga = pd.read_csv('dramaga.csv', sep=',')
data_tang = pd.read_csv('tangerang.csv', sep=',')


data_balai = data_balai.rename(columns={'0':'lat', '0.1':'lon', '0.2':'jam', '0.3':'jenis',})
data_dramaga = data_dramaga.rename(columns={'0':'lat', '0.1':'lon', '0.2':'jam', '0.3':'jenis',})
data_tang = data_tang.rename(columns={'0':'lat', '0.1':'lon', '0.2':'jam', '0.3':'jenis',})

data_balai = data_balai[(data_balai.loc[:,'jam'].str[0:8]=='12:00 AM')]
data_dramaga = data_dramaga[(data_dramaga.loc[:,'jam'].str[0:8]=='12:00 AM')]
data_tang = data_tang[(data_tang.loc[:,'jam'].str[0:8]=='12:00 AM')]

data_balai = data_balai.loc[:,'lat':'lon']
data_dramaga = data_dramaga.loc[:,'lat':'lon']
data_tang = data_tang.loc[:,'lat':'lon']

data_balai = data_balai.values
data_dramaga = data_dramaga.values
data_tang = data_tang.values

X = np.vstack((data_balai, data_dramaga, data_tang))
X = X[:, [1,0]]
clust = OPTICS(min_samples=50, xi=.05, min_cluster_size=.05)

clust.fit(X)

space = np.arange(len(X))
reachability = clust.reachability_[clust.ordering_]
labels = clust.labels_[clust.ordering_]

plt.figure(figsize=(10, 7))
G = gridspec.GridSpec(1, 1)
ax1 = plt.subplot(G[:, :])

# OPTICS
colors = ['g.', 'r.', 'b.', 'y.', 'c.']
for klass, color in zip(range(0, 5), colors):
    Xk = X[clust.labels_ == klass]
    lonrat = Xk[:,0].mean()
    latrat = Xk[:,1].mean()
    print(lonrat)
    print(latrat)
    ax1.plot(lonrat, latrat, color, alpha=0.3)
ax1.plot(X[clust.labels_ == -1, 0], X[clust.labels_ == -1, 1], 'k+', alpha=0.1)
ax1.set_ylim([-10, 0])
ax1.set_xlim([100, 110])
ax1.set_title('Automatic Clustering\nOPTICS')

plt.tight_layout()
plt.show()
