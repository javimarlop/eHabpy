#!/usr/bin/python

"""
Adpated from: https://gist.github.com/geberl/c65517bf8273552486f9a8954e80ddf4

Installation: Get the *Miniconda* Python Distribution - not the Python distrubution from python.org!
- https://conda.io/miniconda.html

Then install modules:
- `cd ~/miniconda3/bin`
- `./conda install numpy pandas matplotlib`

Original source:
- https://stackoverflow.com/questions/24659005/radar-chart-with-multiple-scales-on-multiple-axes
- That code has problems with 5+ axes though
"""

import numpy as np
import matplotlib.pyplot as plt

# Optionally use different styles for the graph
# Gallery: http://tonysyu.github.io/raw_content/matplotlib-style-gallery/gallery.html
# import matplotlib
# matplotlib.style.use('dark_background')  # interesting: 'bmh' / 'ggplot' / 'dark_background'


class Radar(object):
    def __init__(self, figure, title, labels, rect=None):
        if rect is None:
            rect = [0.05, 0.05, 0.9, 0.9]

        self.n = len(title)
        self.angles = np.arange(0, 360, 360.0/self.n)

        self.axes = [figure.add_axes(rect, projection='polar', label='axes%d' % i) for i in range(self.n)]

        self.ax = self.axes[0]
        self.ax.set_thetagrids(self.angles, labels=title, fontsize=12)

        for ax in self.axes[1:]:
            ax.patch.set_visible(False)
            ax.grid(False)
            ax.xaxis.set_visible(False)

        for ax, angle, label in zip(self.axes, self.angles, labels):
            ax.set_rgrids(range(1, 6), angle=angle, labels=label)
            ax.spines['polar'].set_visible(False)
            ax.set_ylim(0, 5.5)

    def plot(self, values, *args, **kw):
        angle = np.deg2rad(np.r_[self.angles, self.angles[0]])
        values = np.r_[values, values[0]]
        self.ax.plot(angle, values, *args, **kw)


if __name__ == '__main__':
    fig = plt.figure(figsize=(6, 6))

    tit = ['THDI','TotEco','SIH','THD','THR']

    lab = [
        ['15','30','45','60','75'],
        list('12345'),
        ['0.3','0.6','0.9','1.2','1.5'],
        ['1.5', '3.0', '4.5', '6.0', '7.5'],
        ['0.15','0.30','0.45','0.60','0.75']
    ]

    df = np.genfromtxt('thdi_res.csv',delimiter=' ', names=True, usecols=(3,4,5,6,7,8))
  
# Kakadu
radar2 = Radar(fig, tit, lab)
n = 5
radar2.plot([df[n,][5]*5/75, df[n,][0], df[n,][2]*5/1.5, df[n,][3]*5/7.5, df[n,][4]*5/0.75],  '-', lw=2, color='b', alpha=0.4, label='first')

fig.suptitle('Canaima', size=16, y = 1, ha = 'right')

fig.savefig('rdplot_can.png', bbox_inches='tight')
