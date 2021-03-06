from matplotlib import pyplot
from shapely.geometry import Polygon
#from descartes.patch import PolygonPatch

#from figures import SIZE

COLOR = {
    True:  '#6699cc',
    False: '#ff3333'
    }

def v_color(ob):
    return COLOR[ob.is_valid]

def plot_coords(ax, ob):
    x, y = ob.xy
    plot(x, y, 'o', color='#999999', zorder=1)
    
#fig = pyplot.figure(1, figsize=SIZE, dpi=90)

# 1: valid polygon
#ax = fig.add_subplot(121)

ext = [(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)]

polygon = Polygon(ext)

pyplot.plot(ext[:,0],ext[:,1])

#plot_coords(ax, polygon.interiors[0])
#plot_coords(ax, polygon.exterior)


pyplot.show()

