from numpy import arange
from scenariogeneration import xodr
from shapely.geometry import LineString, Point
import os
import math
from scipy.interpolate import interp1d, CubicSpline
from scipy import stats, mean

def resample_linear(pts, d=1.0):
    ref = LineString(pts)
    new_pts = []
    n = math.ceil(ref.length / d)
    for idx in range(0, n+1):
        pt = ref.interpolate(idx/n, normalized=True)
        new_pts.append((pt.x, pt.y))
    return new_pts

def resample_cubic(pts, a, b, d):
    param = list(range(0, len(pts)))
    if 1:
        ref_x = CubicSpline(param, [x for x,y in pts], bc_type=((1, math.cos(a)), (1, math.cos(b))), extrapolate=True)
        ref_y = CubicSpline(param, [y for x,y in pts], bc_type=((1, math.sin(a)), (1, math.sin(b))), extrapolate=True)
    else:
        ref_x = CubicSpline(param, [x for x,y in pts], bc_type=((2, 0.0), (2, 0.0)), extrapolate=True)
        ref_y = CubicSpline(param, [y for x,y in pts], bc_type=((2, 0.0), (2, 0.0)), extrapolate=True)
    new_pts = []
    n = math.ceil(param[-1] / d)
    for idx in range(0, n+1):
        x = ref_x(param[-1] * idx / n)
        y = ref_y(param[-1] * idx / n)
        new_pts.append((x, y))
    return new_pts

def pts_dir(pts):
    pts_x = [x for x,y in pts]
    pts_y = [y for x,y in pts]
    min_pts_x = min(pts_x)
    min_pts_y = min(pts_y)
    max_pts_x = max(pts_x)
    max_pts_y = max(pts_y)
    if max_pts_y - min_pts_y < max_pts_x - min_pts_x:
        result = stats.linregress(pts_x, pts_y)
        theta = math.atan(result.slope)
        print("\t heading: %.1f degree" % (theta/math.pi*180))
        theta -= math.pi/2
        return math.cos(theta), math.sin(theta)
    else:
        result = stats.linregress(pts_y, pts_x)
        theta = math.atan(result.slope)
        print("\t heading: %.1f degree" % (theta/math.pi*180))
        theta -= math.pi/2
        return math.cos(theta), math.sin(theta)

def export_road(odr, road, road_id):
    planview = xodr.PlanView()
    s = 0.0
    ref_line = road.ref_line
    pts = [(ref_line[0][idx], ref_line[1][idx]) for idx in range(len(ref_line[0]))]

    ref = LineString(pts)
    #CubicSpline([x for x,y in pts], [y for x,y in pts], bc_type=((2, 0.0), (2, 0.0)))
    for idx in range(len(pts)-1):
        dx = pts[idx+1][0] - pts[idx][0]
        dy = pts[idx+1][1] - pts[idx][1]
        heading = math.atan2(dy, dx)
        dist = math.sqrt(dx*dx + dy*dy)
        planview.add_fixed_geometry(xodr.Line(dist), pts[idx][0], pts[idx][1], heading, s)
        s += dist

    centerlane = xodr.Lane(lane_type=xodr.LaneType.median)
    lanesection = xodr.LaneSection(0, centerlane)
    left_bnd_st = None
    for idx, (lane_subid, lane) in enumerate(road.lanes.items()):
        #if idx > 1:
        #    continue
        width_a = []
        width_b = []
        soffset = []
        right_bnd_pts = [(lane.right_bnd[0][idx], lane.right_bnd[1][idx]) for idx in range(len(lane.right_bnd[0]))]
        # if lane.full_id == "557024172,0,0,36,0":
        #     new_right_bnd_pts = resample_linear(right_bnd_pts, 0.1)
        #     print(f"[{lane.full_id}] end   ","    right_bnd_pts: %.10f %.10f" % (right_bnd_pts[-1][0], right_bnd_pts[-1][1]))
        #     print(f"[{lane.full_id}] end   ","new_right_bnd_pts: %.10f %.10f" % (new_right_bnd_pts[-1][0], new_right_bnd_pts[-1][1]))
        # if lane.full_id == "557024172,0,0,37,0":
        #     new_right_bnd_pts = resample_linear(right_bnd_pts, 0.1)
        #     print(f"[{lane.full_id}] start ","    right_bnd_pts: %.10f %.10f" % (right_bnd_pts[0][0], right_bnd_pts[0][1]))
        #     print(f"[{lane.full_id}] start ","new_right_bnd_pts: %.10f %.10f" % (new_right_bnd_pts[0][0], new_right_bnd_pts[0][1]))

        ## todo: resampling introduces more accumulated errors on lanes of curved roads.
        #right_bnd_pts = resample_linear(right_bnd_pts, 0.1)
        
        #print(lane_subid, right_bnd.length)
        right_bnd_s = []
        right_bnd_t = []
        old_pt2 = None
        for pt_idx, pt in enumerate(right_bnd_pts):
            pt2 = Point(pt[0], pt[1])
            if old_pt2 is not None and old_pt2.distance(pt2) < 1.0:
                continue
            d = ref.distance(pt2)
            s = ref.project(pt2)
            right_bnd_s.append(s)
            right_bnd_t.append(d)
            if left_bnd_st is not None:
                d -= left_bnd_st(s)
            if pt_idx == 0:
                s = 0.0
            elif pt_idx == len(right_bnd_pts) - 1:
                s = LineString(right_bnd_pts).length
            width_a.append(d)
            soffset.append(s)
            old_pt2 = pt2
        for idx in range(len(width_a)-1):
            width_b.append((width_a[idx+1]-width_a[idx])/(soffset[idx+1]-soffset[idx]+0.000001))
        width_b.append(0.0)
        left_bnd_st = interp1d(right_bnd_s, right_bnd_t, fill_value="extrapolate")

        # if lane.full_id == "557024172,0,0,36,0":
        #     print(f"[{lane.full_id}] end   ","    width: %.10f" % (width_a[-1]))
        # if lane.full_id == "557024172,0,0,37,0":
        #     print(f"[{lane.full_id}] start ","    width: %.10f" % (width_a[0]))

        lanesection.add_right_lane(xodr.Lane(lane_type=xodr.LaneType.driving, a=width_a, b=width_b, soffset=soffset))

    # lanesection.add_left_lane(xodr.Lane(lane_type=xodr.LaneType.median, a=[0.3, 2.3, 2.3, 0.3], b=[0.2, 0.0, -0.2, 0.0], soffset=[0.0, 10.0, 20.0, 30.0]))
    # lanesection.add_right_lane(xodr.Lane(lane_type=xodr.LaneType.median, a=0.3))

    # left_lane_with_roadmark = xodr.Lane(a=[4, 5, 5, 4], b=[0.1, 0.0, -0.1, 0.0], soffset=[0.0, 10.0, 20.0, 30.0])
    # left_lane_with_roadmark.add_roadmark(xodr.STD_ROADMARK_BROKEN)
    # lanesection.add_left_lane(left_lane_with_roadmark)
    # right_lane_with_roadmark = xodr.Lane(a=[4, 5, 5, 4], b=[0.1, 0.0, -0.1, 0.0], soffset=[0.0, 10.0, 20.0, 30.0])
    # right_lane_with_roadmark.add_roadmark(xodr.STD_ROADMARK_SOLID)
    # lanesection.add_right_lane(right_lane_with_roadmark)

    lanes = xodr.Lanes()
    lanes.add_lanesection(lanesection)

    road = xodr.Road(road_id, planview, lanes)
    odr.add_road(road)

def export(my_map):
    # create the opendrive
    odr = xodr.OpenDrive("myroad")
    for idx, (road_id, road) in enumerate(my_map.roads.items()):
        #if road_id != "557024172,0,0,66":
        #    continue
        print(road_id)
        export_road(odr, road, idx)
    #odr.adjust_roads_and_lanes()
    odr.write_xml("test.xodr")
    print("Done")