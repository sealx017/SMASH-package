import os

RCFGL_path = '/Users/sealso/Documents/GitHub/SMASH/'
os.chdir(RCFGL_path)

from SMASH import SMASH
import pickle
import seaborn as sns
import pandas as pd


with open("/Users/sealso/Documents/GitHub/SMASH/Data/Mfish.pickle", "rb") as fp:  
 df = pickle.load(fp) 

Gene = df[0]
Cords = df[1]


SMASH_Result = SMASH(Gene, Cords)
SMASH_Result['SMASH'].values[:,1]
SMASH_Result['SPARK-X']


p = Gene .shape[1]

Cords = Cords.set_axis(["X", "Y"], axis = 1)
Data = pd.concat([Cords, Gene['Npy1r']], axis = 1)
Gene_name = 'Npy1r'
pal = sns.color_palette("crest")
pal.as_hex()[3]

sns.color_palette("crest", as_cmap=True)
out = sns.scatterplot(data=Data, x="X", y="Y", hue = Gene_name, s = 3, palette = 'crest')
sns.move_legend(
    out, "lower center",
    bbox_to_anchor=(0.5, 1), ncol=6, title='Relative Expression of ' + Gene_name, frameon=False,
)
out.invert_yaxis()


import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def pattern_plot_SMASH(Data, Gene_names, title_size = 2, s = 0.5, cmap = 'crest'):
    fig, ax = plt.subplots()
    ax.scatter(Data['X'], Data['Y'], s = 0.5, alpha = 1, c = Data[Gene_name], cmap = cmap)
    ax.set_title(label = Gene_name, loc = 'center')
    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size='5%', pad=0.6, pack_start = True)
    fig.add_axes(cax)
    fig.colorbar(ax.imshow, cax = cax, orientation = 'horizontal')
    cbar.set_label('Relative Expression')
    
    
    return ax


plt.scatter(Data['X'], Data['Y'], s = 0.5, alpha = 1, c = Data[Gene_name], cmap = 'crest')
cbar = plt.colorbar()



def Expression_plot(Data, Gene_name, s = 0.5, cmap = 'crest'):
    plt.figure()
    ax = plt.gca()
    ax.set_title(label = Gene_name, loc = 'center')
    im = ax.scatter(Data['X'], Data['Y'], s = 0.5, alpha = 1, c = Data[Gene_name]/max(Data[Gene_name]), cmap = cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)     
    plt.colorbar(im, cax=cax)
    return None


Expression_plot(Data, Gene_name, s = 0.5, cmap = 'crest')
