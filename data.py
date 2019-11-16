import pandas as pd
import matplotlib.pyplot as plt

data_balai = pd.read_csv('balai.csv', sep=',')
data_dramaga = pd.read_csv('dramaga.csv', sep=',')
data_tang	= pd.read_csv('tangerang.csv', sep=',')

data_balai = data_balai.rename(columns={'0':'lat', '0.1':'lon', '0.2':'jam', '0.3':'jenis',})
data_dramaga = data_dramaga.rename(columns={'0':'lat', '0.1':'lon', '0.2':'jam', '0.3':'jenis',})
data_tang = data_tang.rename(columns={'0':'lat', '0.1':'lon', '0.2':'jam', '0.3':'jenis',})

data_balai = data_balai[(data_balai.loc[:,'jam'].str[0:8]=='12:30 AM')]
data_dramaga = data_dramaga[(data_dramaga.loc[:,'jam'].str[0:8]=='12:30 AM')]
data_tang = data_tang[(data_tang.loc[:,'jam'].str[0:8]=='12:30 AM')]

#dramaga 6°33'11.3"S 106°44'34.8"E
#balai II 6°18'09.2"S 106°45'22.9"E
#tan 6°10'16.7"S 106°38'48.2"E

balai = (data_balai['lat'], data_balai['lon'])
dramaga = (data_dramaga['lat'], data_dramaga['lon'])
tang = (data_tang['lat'], data_tang['lon'])

gab = (balai, dramaga, tang)
colors = ('red', 'green', 'blue')
groups = ('balai', 'dramaga', 'tangerang')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for data, color, group in zip(gab, colors, groups):
    x,y = data
    ax.scatter(y,x, alpha=0.5, c=color, edgecolors='none', s=10, label=group)

plt.ylim(-10,0)
plt.xlim(100,120)
plt.title('20190101  12:00 AM')
plt.legend(loc=2)
plt.show()
