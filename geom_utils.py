import math
import matplotlib.pyplot as plt

class WorldBox:
    x_min, x_max, y_min, y_max = 1e30, -1e30, 1e30, -1e30
    @classmethod
    def update(cls, xx, yy):
        cls.x_min = min(cls.x_min, min(xx))
        cls.x_max = max(cls.x_max, max(xx))
        cls.y_min = min(cls.y_min, min(yy))
        cls.y_max = max(cls.y_max, max(yy))

    @classmethod
    def get(cls):
        return cls.x_min, cls.x_max, cls.y_min, cls.y_max

    @classmethod
    def update_fig_range(cls, ax, w, h):
        x_min, x_max, y_min, y_max = cls.get()
        dx = x_max - x_min
        dy = y_max - y_min
        if dx/dy > w/h:
            dy = dx/w*h
            y = (y_max + y_min) / 2
            y_min = y - dy / 2
            y_max = y + dy / 2
        else:
            dx = dy*w/h
            x = (x_min + x_max) / 2
            x_min = x - dx / 2
            x_max = x + dx / 2
        ax.set_xlim((x_min, x_max))
        ax.set_ylim((y_min, y_max))        
        mng = plt.get_current_fig_manager()
        mng.resize(w, h)

def first(c):
  return next(iter(c))

def last(c):
  return next(iter(reversed(c)))

def dist2(p1, p2):
  return ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

def is_almost_the_same_pt(p1, p2):
  return dist2(p1, p2) < 0.05 * 0.05

def curvature(x, y):
  triangle_area = abs((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0]))
  a = dist2((x[0], y[0]), (x[1], y[1]))
  b = dist2((x[2], y[2]), (x[1], y[1]))
  c = dist2((x[0], y[0]), (x[2], y[2]))
  return 4 * triangle_area / math.sqrt(a * b * c)

def pt_hash(x, y):
  return "%.1f,%.1f" % (x, y)

def is_almost_the_same_pt2(p1, p2):
  return pt_hash(*p1) == pt_hash(*p2)

def xxyy2xyxy(xxyy):
  return [(xxyy[0][idx], xxyy[1][idx]) for idx in range(len(xxyy[0]))]

def xyxy2xxyy(xy):
  return ([x for x, y in xy], [y for x, y in xy])