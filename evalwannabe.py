from sklearn.cluster import OPTICS, cluster_optics_dbscan
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# DATA NEXSTRORM
data_nexstorm = pd.read_csv('balai.csv', sep=',')
data_nexstorm = data_balai[(data_balai.loc[:,'datetime_utc'].str[11:13]=='15')&(data_balai.loc[:,'datetime_utc'].str[14:16]>'0')&(data_balai.loc[:,'datetime_utc'].str[14:16]<'15')]
data_nexstorm = data_balai.rename(columns={'latitude':'lat', 'longitude':'lon'})
data_nexstorm['sensor'] = 'namasensor'
data_nexstorm = data_balai.loc[:,['lat','lon','sensor']]

#DATA LD2000
data_balai = pd.read_csv('balai.csv', sep=',')
data_dramaga = pd.read_csv('dramaga.csv', sep=',')
data_tang = pd.read_csv('tangerang.csv', sep=',')

data_balai = data_balai.rename(columns={'0':'lat', '0.1':'lon', '0.2':'jam', '0.3':'jenis'})
data_dramaga = data_dramaga.rename(columns={'0':'lat', '0.1':'lon', '0.2':'jam', '0.3':'jenis'})
data_tang = data_tang.rename(columns={'0':'lat', '0.1':'lon', '0.2':'jam', '0.3':'jenis'})

data_balai['sensor']='Balai'
data_dramaga['sensor']='Dramaga'
data_tang['sensor']='Tangerang'

data_balai = data_balai[(data_balai.loc[:,'jam'].str[0:8]=='12:00 AM')]
data_dramaga = data_dramaga[(data_dramaga.loc[:,'jam'].str[0:8]=='12:00 AM')]
data_tang = data_tang[(data_tang.loc[:,'jam'].str[0:8]=='12:00 AM')]

data_balai = data_balai.loc[:,['lat','lon','sensor']]
data_dramaga = data_dramaga.loc[:,['lat','lon','sensor']]
data_tang = data_tang.loc[:,['lat','lon','sensor']]

kabeh_data = pd.merge(data_balai, data_dramaga, how='outer')
kabeh_data = pd.merge(kabeh_data, data_tang, how='outer')
kabeh_data = pd.merge(kabeh_data, data_nexstorm, how='outer')

data_balai = data_balai.values
data_dramaga = data_dramaga.values
data_tang = data_tang.values

X = np.vstack((data_balai, data_dramaga, data_tang))
X = X[:, [1,0]]
clust = OPTICS(min_samples=50, xi=.05, min_cluster_size=.05)

clust.fit(X[:, [1,0]])

space = np.arange(len(X))
reachability = clust.reachability_[clust.ordering_]
labels = clust.labels_[clust.ordering_]

plt.figure(figsize=(10, 7))
G = gridspec.GridSpec(2, 3)
ax1 = plt.subplot(G[0, :])
ax2 = plt.subplot(G[1, :])

# Reachability plot
#colors = ['g.', 'r.', 'b.']
colors = ['g.', 'r.', 'b.', 'y.', 'c.']
for klass, color in zip(range(0, 5), colors):
    Xk = space[labels == klass]
    Rk = reachability[labels == klass]
    ax1.plot(Xk, Rk, color, alpha=0.3)
ax1.set_ylim([0, 4])
ax1.plot(space[labels == -1], reachability[labels == -1], 'k.', alpha=0.3)
ax1.plot(space, np.full_like(space, 2., dtype=float), 'k-', alpha=0.5)
ax1.plot(space, np.full_like(space, 0.5, dtype=float), 'k-.', alpha=0.5)
ax1.set_ylabel('Reachability (epsilon distance)')
ax1.set_title('20190101 12:00 AM')

# OPTICS
colors = ['g.', 'r.', 'b.', 'y.', 'c.']

for klass, color in zip(range(0, 5), colors):
    Xk = X[clust.labels_ == klass]
    kepilih = pd.DataFrame(Xk)
    kepilih = kepilih.rename(columns={0:'lon', 1:'lat'})
    kepilih = kepilih.astype('float64')
    kepilih_final = pd.merge(kepilih, kabeh_data, how='outer', on='lon')
    kepilih_final = kepilih_final.dropna()
    kepilih_final = kepilih_final['sensor'].unique().tolist()
    print(kepilih_final)
    lonrat = Xk[:,0].mean()
    latrat = Xk[:,1].mean()
    print(lonrat)
    print(latrat)
    ax2.plot(lonrat, latrat, color, alpha=0.3)
ax2.plot(X[clust.labels_ == -1, 0], X[clust.labels_ == -1, 1], 'k+', alpha=0.1)
ax2.set_ylim([-10, 0])
ax2.set_xlim([100, 110])
ax2.set_title('Automatic Clustering\nOPTICS')

plt.tight_layout()
plt.show()
