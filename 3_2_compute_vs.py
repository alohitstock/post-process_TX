import matplotlib.pyplot as plt 

def motion(event):  
    x = event.xdata
    y = event.ydata
    if event.inaxes == axes[0]:
        ln_1.set_data(x,y)
    if event.inaxes == axes[1]:
        ln_2.set_data(x,y)
    plt.draw()

_, axes = plt.subplots(2, 1)
ln_1, = axes[0].plot([],[],'o')

ln_2, = axes[1].plot([],[],'x')
print(ln_2)

plt.connect('motion_notify_event', motion)
plt.show()