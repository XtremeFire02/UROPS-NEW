import os, math, json
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np, pandas as pd, matplotlib.pyplot as plt

from matplotlib import font_manager

wayPtFile = "wayPtSantosEnd.csv"
coastColor = "gray"
pathColor = "dimgray"

if (os.path.isfile(wayPtFile) == False):
    print(f"File {wayPtFile} not found")
    exit()

#def shiftLon360(lon):
#    lonS = np.where(lon > 180, lon - 360, lon)
#    return lonS

#wObj = np.load(windUVFile)

wayPtDF = pd.read_csv(wayPtFile)

wayPtLat = np.array(wayPtDF["lat"])
wayPtLon = np.array(wayPtDF["lon"])


#proj = ccrs.Gnomonic()
proj = ccrs.PlateCarree()

fig = plt.figure(figsize=(10, 5))
ax = plt.axes(projection=proj)

#ax.set_extent([lonLb, lonUb, latLb, latUb], crs=proj)

ax.coastlines(resolution="50m", color=coastColor)

#ly = np.arange(math.ceil(latLb), math.floor(latUb) + 1, 4)
#lx = np.arange(math.ceil(lonLb), math.floor(lonUb) + 1, 4)
#xlocs=lx, ylocs=ly,

gl = ax.gridlines(crs=proj, draw_labels=True, linewidth=1, color=coastColor, alpha=0.5, linestyle='--')

gl.top_labels = False
gl.right_labels = False

plt.title("File: " + wayPtFile)

plt.plot(wayPtLon, wayPtLat, c=pathColor, lw=2, linestyle=":")
plt.annotate("start", xy=(wayPtLon[0], wayPtLat[0]))
plt.annotate("end", xy=(wayPtLon[-1], wayPtLat[-1]))

plt.show()
