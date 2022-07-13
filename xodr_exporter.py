from collections import OrderedDict
from numpy import arange
from scenariogeneration import xodr
from shapely.geometry import LineString, Point
import os
import math
from scipy.interpolate import interp1d, CubicSpline
from scipy import stats, mean
from geom_utils import xxyy2xyxy, xyxy2xxyy


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
    xx, yy = xyxy2xxyy(pts)
    if 1:
        ref_x = CubicSpline(param, xx, bc_type=((1, math.cos(a)), (1, math.cos(b))), extrapolate=True)
        ref_y = CubicSpline(param, yy, bc_type=((1, math.sin(a)), (1, math.sin(b))), extrapolate=True)
    else:
        ref_x = CubicSpline(param, xx, bc_type=((2, 0.0), (2, 0.0)), extrapolate=True)
        ref_y = CubicSpline(param, yy, bc_type=((2, 0.0), (2, 0.0)), extrapolate=True)
    new_pts = []
    n = math.ceil(param[-1] / d)
    for idx in range(0, n+1):
        x = ref_x(param[-1] * idx / n)
        y = ref_y(param[-1] * idx / n)
        new_pts.append((x, y))
    return new_pts

def pts_dir(pts):
    pts_x, pts_y = xyxy2xxyy(pts)
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

def export_road(odr, road):
    planview = xodr.PlanView()
    s = 0.0
    pts = xxyy2xyxy(road.ref_line)

    ref = LineString(pts)
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
    lane_offset_s = []
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

        ## TODO: Resampling introduces more accumulated errors on lanes of curved roads.
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
        if lane.is_fake():
            lane_offset_s = soffset
            lane_offset_a = width_a
            lane_offset_b = width_b
            continue

        new_xodr_lane = xodr.Lane(lane_type=xodr.LaneType.driving, a=width_a, b=width_b, soffset=soffset)
        if len(lanesection.rightlanes) > 0:
            roadmark_dashed = xodr.RoadMark(xodr.RoadMarkType.broken, 0.01)
            lanesection.rightlanes[-1].add_roadmark(roadmark_dashed)
        lanesection.add_right_lane(new_xodr_lane)
        lane.xodr = new_xodr_lane

    lanes = xodr.Lanes()
    if len(lane_offset_s) > 0:
        for s, a, b in zip(lane_offset_s, lane_offset_a, lane_offset_b):
            lanes.add_laneoffset(xodr.LaneOffset(s, -a, -b))
    lanes.add_lanesection(lanesection)

    junction_id = road.junction.id if road.junction is not None else -1
    xodr_road = xodr.Road(road.short_id, planview, lanes, road_type=junction_id)
    xodr_road.original_id = road.id
    odr.add_road(xodr_road)
    road.xodr = xodr_road

def export_road_linkage(odr, road_a):
    road_set = set()

    for idx in range(2):
        for linkage in road_a.linkage[idx]:
            road_b, start_or_end_b = linkage
            road_set.add(road_b)
            if start_or_end_b == "end":
                #road_b.xodr.add_successor(xodr.ElementType.road, road_a.xodr.id, xodr.ContactPoint.start)
                road_a.xodr.add_predecessor(xodr.ElementType.road, road_b.xodr.id, xodr.ContactPoint.end)
            elif start_or_end_b == "start":
                #road_b.xodr.add_predecessor(xodr.ElementType.road, road_a.xodr.id, xodr.ContactPoint.end)
                road_a.xodr.add_successor(xodr.ElementType.road, road_b.xodr.id, xodr.ContactPoint.start)
            else: # "junction"
                junction_id = road_b
                if idx == 0:
                    road_a.xodr.add_predecessor(xodr.ElementType.junction, junction_id)
                else:
                    road_a.xodr.add_successor(xodr.ElementType.junction, junction_id)

    for lane_id, lane in road_a.lanes.items():
        for predecessor in lane.predecessors:
            for linkage in road_a.linkage[0]:
                from_junction = linkage is not None and linkage[1] == "junction"
                if from_junction or predecessor.road in road_set:
                    lane.xodr.add_link("predecessor", predecessor.xodr.lane_id)
        for successor in lane.successors:
            for linkage in road_a.linkage[1]:
                to_junction = linkage is not None and linkage[1] == "junction"
                if to_junction or successor.road in road_set:
                    lane.xodr.add_link("successor", successor.xodr.lane_id)

def export_default_junction(odr, junction):
    j = xodr.Junction(junction.set_id(), junction.id)
    connections = OrderedDict()
    for road in junction.connecting_roads:
        for lane in road.lanes.values():
            for prev_lane in lane.predecessors:
                connection_id = (prev_lane.road.short_id, lane.road.short_id)
                if connection_id not in connections:
                    connections[connection_id] = []
                connections[connection_id].append((prev_lane.xodr.lane_id, lane.xodr.lane_id))
    for (incoming_road_id, connecting_road_id), links in connections.items():
        c = xodr.Connection(incoming_road_id, connecting_road_id, xodr.ContactPoint.start)
        for incoming_lane_id, connecting_lane_id in links:
            c.add_lanelink(incoming_lane_id, connecting_lane_id)
        j.add_connection(c)
    odr.add_junction(j)

def export(xodr_filename, my_map, offset_x, offset_y, map_ver, georef):
    odr = xodr.OpenDrive(xodr_filename.replace(".json", ""))
    #odr.add_offset(offset_x, offset_y, 0.0)
    #odr.add_user_data("inceptioMapVersion", map_ver)
    #odr.add_user_data("inceptioMapGeoRef", georef)

    for road_id, road in my_map.roads.items():
        print(f"Exporting Road[{road_id}]...")
        export_road(odr, road)
        
    for road_id, road in my_map.roads.items():
        print(f"Exporting Linkages of Road[{road_id}]...")
        export_road_linkage(odr, road)

    for junction in my_map.default_junctions:
        print(f"Exporting Default Junction[{junction.set_id()}]...")
        export_default_junction(odr, junction)

    print("Write the OpenDRIVE file...")
    odr.write_xml(xodr_filename)
    print("Done")