import math
import matplotlib.pyplot as plt

class WorldBox:
    x_min, x_max, y_min, y_max = 1e30, -1e30, 1e30, -1e30

    @classmethod
    def clear(cls):
      cls.x_min, cls.x_max, cls.y_min, cls.y_max = 1e30, -1e30, 1e30, -1e30

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

def curvature(pts):
  triangle_area = abs((pts[1][0]-pts[0][0])*(pts[2][1]-pts[0][1]) - (pts[1][1]-pts[0][1])*(pts[2][0]-pts[0][0]))
  a = dist2(pts[0], pts[1])
  b = dist2(pts[2], pts[1])
  c = dist2(pts[0], pts[2])
  return 4 * triangle_area / math.sqrt(a * b * c)

def pt_hash(pt):
  return "%.1f,%.1f" % (pt[0], pt[1])

def is_almost_the_same_pt2(p1, p2):
  return pt_hash(*p1) == pt_hash(*p2)

def xxyy2xyxy(xxyy):
  return [(xxyy[0][idx], xxyy[1][idx]) for idx in range(len(xxyy[0]))]

def xyxy2xxyy(xy):
  return ([x for x, y in xy], [y for x, y in xy])

def clip_xyxy(xyxy, length):
  if not isinstance(xyxy, list):
    xyxy = list(xyxy)
  assert(len(xyxy) > 1)
  d = 0.0
  xyxy2 = [xyxy[0]]
  for idx in range(1, len(xyxy)):
    xy = xyxy[idx]
    dd = math.sqrt(dist2(xyxy2[-1], xy))
    xyxy2.append(xy)
    d += dd
    if d >= length:
      break
  return xyxy2
  