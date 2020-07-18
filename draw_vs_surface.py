import plotly.graph_objects as go
import numpy as np
import plotly.offline as po

# prepare Vs Surface
data_x = np.linspace(0, 100, 1000)
data_y = np.linspace(-50, 50, 2001)
data_surface1 = np.zeros((1500, 1001))
for x in range(data_surface1.shape[0]):
    for y in range(data_surface1.shape[1]):
        if (2*y - (x - 500)) < 0 or (2*y + (x - 500)) < 0:
            data_surface1[x, y] = 0
        else:
            data_surface1[x, y] =(((2*y + (x - 500)) / 2)**(1/3) * ((2*y - (x - 500)) / 2)**(2/3)) ** (0.7)

fig = go.Figure(data=[go.Surface(x=data_x, y=data_y, z=data_surface1)])

fig.update_layout(title='Vs Surface', autosize=False,
                  width=900, height=900,
                  margin=dict(l=65, r=50, b=65, t=90))

fig.show()
