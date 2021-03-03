import sympy
from sympy.abc import x, y
import numpy as np
import matplotlib.pyplot as plt


def cylinder_stream_function():
    r = sympy.sqrt(x**2 + y**2)
    theta = sympy.atan2(y, x)
    return -2*r*sympy.cos(theta)-theta-1/(4*sympy.pi)*sympy.log(r)

def velocity_field(psi):
    u = sympy.lambdify((x, y), psi.diff(y), 'numpy')
    v = sympy.lambdify((x, y), -psi.diff(x), 'numpy')
    return u, v

def plot_streamlines(ax, u, v, xlim=(-3, 3), ylim=(-3, 3)):
    x0, x1 = xlim
    y0, y1 = ylim
    Y, X =  np.ogrid[y0:y1:10j, x0:x1:10j] 
    stream_points = np.array(list(zip(np.arange(-3,3,.2), -3*np.ones(30))))
    ax.streamplot(X, Y, u(X, Y), v(X, Y),color='blue',start_points=stream_points,density=30)
    
def format_axes(ax):
    ax.set_aspect('equal')
    ax.figure.subplots_adjust(bottom=0, top=1, left=0, right=1)

    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.axis('off')
        
psi = cylinder_stream_function()
u, v = velocity_field(psi)

xlim = ylim = (-3, 3)
fig, ax = plt.subplots(figsize=(4, 4))
plot_streamlines(ax, u, v, xlim, ylim,)

format_axes(ax)
fig.savefig('straumlinur.eps', format='eps', dpi=1200)