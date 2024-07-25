#using simplify to polygon segments
#using earcut to accelarate make edges
#filiter vertex > 5000
#using skel-algorithm --0723
import math
from collections import defaultdict
import heapq
import glue
import export_to_drc
import os
import simplify
import earcut
import random
import polyskel
def generate_drc_cal(layout_path, layout_primary, drc_results_database, rule_name, file_name):
    content = """LAYOUT PATH  "{}"
//DRC RESULTS DATABASE PRECISION 0.1
LAYOUT PRIMARY "{}"
LAYOUT SYSTEM GDSII
DRC RESULTS DATABASE "{}" ASCII PSEUDO
DRC MAXIMUM RESULTS ALL
DRC MAXIMUM VERTEX ALL
DRC CELL NAME YES CELL SPACE XFORM
VIRTUAL CONNECT COLON YES
VIRTUAL CONNECT REPORT NO
VIRTUAL CONNECT NAME ?
DRC SELECT CHECK
    "M1_polygons"
    "M2_polygons"
    "M3_polygons"
    "M4_polygons"
    "M5_polygons"
    "V1C_polygons"
    "V1P_polygons"
    "V2C_polygons"
    "V2P_polygons"
    "V3C_polygons"
    "V3P_polygons"
    "V4C_polygons"
    "V4P_polygons"
    "Net_polygons"
DRC ICSTATION YES
INCLUDE "/workspace/projs/tah/proj_lib/L13450163/DB_A1/CAD/DRC_LVS/DRC_RUNSET/tah_drc.cal"
INCLUDE "{}"
""".format(layout_path, layout_primary, drc_results_database, rule_name)
    with open(file_name, 'w') as file:
        file.write(content)

def generate_checkbuffer_rule(length_limit, area_limit, perimeter_limit, ground_names, file_name):
    quoted_ground_names = ['"' + name + '"' for name in ground_names]
    ground_names_str = " ".join(quoted_ground_names)
    text = """variable Length_Limit {}
variable Area_Limit {}
variable Perimeter_Limit {}
variable GROUD_NAME\t {}
M2\t\t\t= M2SP OR M2DP
M3\t\t\t= copy M3SP
P2_NET_WITH_GATE\t= P2 INTERACT P2_GATE
PC1_OR_VIA0.COMBINED\t= PC1 OR VIA0.COMBINED
PC1_OR_VIA0.AA\t\t= PC1_OR_VIA0.COMBINED NOT P2
PEPPER_VIA\t\t= ((BS AND BVIA) AND TVIA) AND TS
VIA1.COMBINED_ARRAY\t= V1_ARRAY OR V1C_ARRAY
VIA0.COMBINED_ARRAY\t= OR V0DP V0ACS V0SS V0P_ARRAY V0P_TAC V0P_TSC
PC1.COMBINED_ARRAY\t= OR PC1_ARRAY PC1_TSC PC1_TAB
PC1_VIA0.COMBINED\t= PC1.COMBINED_ARRAY AND VIA0.COMBINED_ARRAY
TSC_PC1_VIA0.COMBINED\t= TSC AND PC1_VIA0.COMBINED
M1_ARRAY_ALL\t\t= M1_ARRAY OR M1DP

//CONNECT P2_GATE P2_NET_WITH_GATE P2
//CONNECT M1 P2_NET_WITH_GATE BY PC1_OR_VIA0.COMBINED
//CONNECT M1 P2 BY PC1_OR_VIA0.COMBINED
//CONNECT M2 M1 BY VIA1.COMBINED
//CONNECT M3SP M2 BY VIA2.COMBINED
//CONNECT M4 M3SP BY VIA3.COMBINED
//CONNECT M5 M4 BY VIA4.COMBINED
//CONNECT TM M5 by TV

// SD connectivity code
CONNECT N_TAP LVNW_CONNECTS BY SD_CT MASK
CONNECT NHV_TAP HVNW_CONNECTS BY SD_CT MASK
CONNECT NHV_TAP HVNW2_CONNECTS BY SD_CT MASK
CONNECT NHV_TAP HVNW_B_CONNECTS BY SD_CT MASK
CONNECT NHV_TAP DNW_CONNECTS BY SD_CT MASK
CONNECT P_TAP LVPW_CONNECTS BY SD_CT MASK
CONNECT PHV_TAP TPW_CONNECTS BY SD_CT MASK
CONNECT SD.COMBINED_NOT_P2 N_TAP NHV_TAP P_TAP PHV_TAP

// P2 connectivity code
CONNECT P2 P2CONNECTS MASK

// M1 connectivity code
TEXT LAYER TEXTM1
TEXT LAYER IOTEXTM1
PORT LAYER TEXT IOTEXTM1
ATTACH TEXTM1 M1CONNECTS MASK
ATTACH IOTEXTM1 M1CONNECTS MASK
//CONNECT M1 M1CONNECTS MASK
CONNECT M1CONNECTS P2CONNECTS BY P2_CT MASK
CONNECT M1CONNECTS SD.COMBINED_NOT_P2 BY SD_CT MASK

// M2 connectivity code
TEXT LAYER TEXTM2
TEXT LAYER IOTEXTM2
PORT LAYER TEXT IOTEXTM2
ATTACH TEXTM2 M2CONNECTS MASK
ATTACH IOTEXTM2 M2CONNECTS MASK
//CONNECT M2SP M2CONNECTS MASK
CONNECT M2CONNECTS M1CONNECTS BY VIA1.COMBINED MASK
CONNECT M2CONNECTS M1CONNECTS BY V1P MASK
CONNECT M2CONNECTS M1CONNECTS BY V1C MASK

// M3 connectivity code
TEXT LAYER TEXTM3
TEXT LAYER IOTEXTM3
PORT LAYER TEXT IOTEXTM3
ATTACH TEXTM3 M3CONNECTS MASK
ATTACH IOTEXTM3 M3CONNECTS MASK
//CONNECT M3SP M3CONNECTS MASK
CONNECT M3CONNECTS M2CONNECTS BY VIA2.COMBINED MASK
CONNECT M3CONNECTS M2CONNECTS BY V2P MASK
CONNECT M3CONNECTS M2CONNECTS BY V2C MASK

// M4 connectivity code
TEXT LAYER TEXTM4
TEXT LAYER IOTEXTM4
PORT LAYER TEXT IOTEXTM4
ATTACH TEXTM4 M4CONNECTS MASK
ATTACH IOTEXTM4 M4CONNECTS MASK
CONNECT M4CONNECTS M3CONNECTS BY VIA3.COMBINED MASK
CONNECT M4CONNECTS M3CONNECTS BY V3P MASK
CONNECT M4CONNECTS M3CONNECTS BY V3C MASK

// M5 connectivity code
TEXT LAYER TEXTM5
TEXT LAYER IOTEXTM5
PORT LAYER TEXT IOTEXTM5
ATTACH TEXTM5 M5CONNECTS MASK
ATTACH IOTEXTM5 M5CONNECTS MASK
CONNECT M5CONNECTS M4CONNECTS BY VIA4.COMBINED MASK
CONNECT M5CONNECTS M4CONNECTS BY V4P MASK
CONNECT M5CONNECTS M4CONNECTS BY V4C MASK

// TM connectivity code
TEXT LAYER TEXTTM
TEXT LAYER IOTEXTTM
PORT LAYER TEXT IOTEXTTM
ATTACH TEXTTM TM MASK
ATTACH IOTEXTTM TM MASK
CONNECT TM M5CONNECTS BY TV MASK

M1_TO_VCC = OR (NET AREA RATIO M1CONNECTS M1_WI_H > 0) (NET AREA RATIO M1CONNECTS M1_WI_L > 0) (NET AREA RATIO M1CONNECTS M2_WI_H > 0) (NET AREA RATIO M1CONNECTS M2_WI_L > 0) (NET AREA RATIO M1CONNECTS M3_WI_H > 0) (NET AREA RATIO M1CONNECTS M3_WI_L > 0) (NET AREA RATIO M1CONNECTS M4_WI_H > 0) (NET AREA RATIO M1CONNECTS M4_WI_L > 0) (NET AREA RATIO M1CONNECTS M5_WI_H > 0) (NET AREA RATIO M1CONNECTS M5_WI_L > 0) (NET AREA RATIO M1CONNECTS TM_WI_H > 0) (NET AREA RATIO M1CONNECTS TM_WI_L > 0)
M1_NOT_TO_VCC = (M1CONNECTS not M1_TO_VCC) not (NET AREA RATIO M1CONNECTS TM > 0)

M2_TO_VCC = OR (NET AREA RATIO M2CONNECTS M1_WI_H > 0) (NET AREA RATIO M2CONNECTS M1_WI_L > 0) (NET AREA RATIO M2CONNECTS M2_WI_H > 0) (NET AREA RATIO M2CONNECTS M2_WI_L > 0) (NET AREA RATIO M2CONNECTS M3_WI_H > 0) (NET AREA RATIO M2CONNECTS M3_WI_L > 0) (NET AREA RATIO M2CONNECTS M4_WI_H > 0) (NET AREA RATIO M2CONNECTS M4_WI_L > 0) (NET AREA RATIO M2CONNECTS M5_WI_H > 0) (NET AREA RATIO M2CONNECTS M5_WI_L > 0) (NET AREA RATIO M2CONNECTS TM_WI_H > 0) (NET AREA RATIO M2CONNECTS TM_WI_L > 0)
M2_NOT_TO_VCC = (M2CONNECTS not M2_TO_VCC) not (NET AREA RATIO M2CONNECTS TM > 0)

M3_TO_VCC = OR (NET AREA RATIO M3CONNECTS M1_WI_H > 0) (NET AREA RATIO M3CONNECTS M1_WI_L > 0) (NET AREA RATIO M3CONNECTS M2_WI_H > 0) (NET AREA RATIO M3CONNECTS M2_WI_L > 0) (NET AREA RATIO M3CONNECTS M3_WI_H > 0) (NET AREA RATIO M3CONNECTS M3_WI_L > 0) (NET AREA RATIO M3CONNECTS M4_WI_H > 0) (NET AREA RATIO M3CONNECTS M4_WI_L > 0) (NET AREA RATIO M3CONNECTS M5_WI_H > 0) (NET AREA RATIO M3CONNECTS M5_WI_L > 0) (NET AREA RATIO M3CONNECTS TM_WI_H > 0) (NET AREA RATIO M3CONNECTS TM_WI_L > 0)
M3_NOT_TO_VCC = (M3CONNECTS not M3_TO_VCC) not (NET AREA RATIO M3CONNECTS TM > 0)

M4_TO_VCC = OR (NET AREA RATIO M4CONNECTS M1_WI_H > 0) (NET AREA RATIO M4CONNECTS M1_WI_L > 0) (NET AREA RATIO M4CONNECTS M2_WI_H > 0) (NET AREA RATIO M4CONNECTS M2_WI_L > 0) (NET AREA RATIO M4CONNECTS M3_WI_H > 0) (NET AREA RATIO M4CONNECTS M3_WI_L > 0) (NET AREA RATIO M4CONNECTS M4_WI_H > 0) (NET AREA RATIO M4CONNECTS M4_WI_L > 0) (NET AREA RATIO M4CONNECTS M5_WI_H > 0) (NET AREA RATIO M4CONNECTS M5_WI_L > 0) (NET AREA RATIO M4CONNECTS TM_WI_H > 0) (NET AREA RATIO M4CONNECTS TM_WI_L > 0)
M4_NOT_TO_VCC = (M4CONNECTS not M4_TO_VCC) not (NET AREA RATIO M4CONNECTS TM > 0)

M5_TO_VCC = OR (NET AREA RATIO M5CONNECTS M1_WI_H > 0) (NET AREA RATIO M5CONNECTS M1_WI_L > 0) (NET AREA RATIO M5CONNECTS M2_WI_H > 0) (NET AREA RATIO M5CONNECTS M2_WI_L > 0) (NET AREA RATIO M5CONNECTS M3_WI_H > 0) (NET AREA RATIO M5CONNECTS M3_WI_L > 0) (NET AREA RATIO M5CONNECTS M4_WI_H > 0) (NET AREA RATIO M5CONNECTS M4_WI_L > 0) (NET AREA RATIO M5CONNECTS M5_WI_H > 0) (NET AREA RATIO M5CONNECTS M5_WI_L > 0) (NET AREA RATIO M5CONNECTS TM_WI_H > 0) (NET AREA RATIO M5CONNECTS TM_WI_L > 0)
M5_NOT_TO_VCC = (M5CONNECTS not M5_TO_VCC) not (NET AREA RATIO M5CONNECTS TM > 0)

area_prop = DFM PROPERTY NET M1_NOT_TO_VCC M2_NOT_TO_VCC M3_NOT_TO_VCC M4_NOT_TO_VCC M5_NOT_TO_VCC
[net_area = area(M1_NOT_TO_VCC) + area(M2_NOT_TO_VCC) + area(M3_NOT_TO_VCC) + area(M4_NOT_TO_VCC) + area(M5_NOT_TO_VCC)]>Area_Limit
[net_peri = perimeter(M1_NOT_TO_VCC) + perimeter(M2_NOT_TO_VCC) + perimeter(M3_NOT_TO_VCC) + perimeter(M4_NOT_TO_VCC) + perimeter(M5_NOT_TO_VCC)]>Perimeter_Limit
P2_A = DFM PROPERTY P2 area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
M1_A = DFM PROPERTY M1_NOT_TO_VCC area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
M2_A = DFM PROPERTY M2_NOT_TO_VCC area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
M3_A = DFM PROPERTY M3_NOT_TO_VCC area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
M4_A = DFM PROPERTY M4_NOT_TO_VCC area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
M5_A = DFM PROPERTY M5_NOT_TO_VCC area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
V1C_A = DFM PROPERTY V1C area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
V1P_A = DFM PROPERTY V1P area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
V2C_A = DFM PROPERTY V2C area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
V2P_A = DFM PROPERTY V2P area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
V3C_A = DFM PROPERTY V3C area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
V3P_A = DFM PROPERTY V3P area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
V4C_A = DFM PROPERTY V4C area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0
V4P_A = DFM PROPERTY V4P area_prop NODAL MULTI [net_area = COUNT(area_prop)]>0

'P2_polygons'{{copy P2_A}}
'M1_polygons'{{copy M1_A}}
'M2_polygons'{{copy M2_A}}
'M3_polygons'{{copy M3_A}}
'M4_polygons'{{copy M4_A}}
'M5_polygons'{{copy M5_A}}
'V1C_polygons'{{copy V1C_A}}
'V1P_polygons'{{copy V1P_A}}
'V2C_polygons'{{copy V2C_A}}
'V2P_polygons'{{copy V2P_A}}
'V3C_polygons'{{copy V3C_A}}
'V3P_polygons'{{copy V3P_A}}
'V4C_polygons'{{copy V4C_A}}
'V4P_polygons'{{copy V4P_A}}
'Net_polygons'{{ornet(ornet(ornet(ornet(ornet P2_A M1_A) M2_A) M3_A) M4_A) M5_A}}
""".format(length_limit, area_limit, perimeter_limit, ground_names_str)
    with open(file_name, 'w') as file:
        file.write(text)


def point_in_polygon(point, polygon):
    x, y = point
    n = len(polygon)
    inside = False

    p1x, p1y = polygon[0]
    for i in range(n):
        p2x, p2y = polygon[(i + 1) % n]
        if ((p1y == y and p2y == y and x >= min(p1x, p2x) and x <= max(p1x, p2x)) or
            (p1x == x and p2x == x and y >= min(p1y, p2y) and y <= max(p1y, p2y))):
            inside = True
            break
        if (p1y == p2y and y == p1y) and ((p1x <= x <= p2x) or (p2x <= x <= p1x)):
            inside = True
            break
        if ((y > min(p1y, p2y)) and (y <= max(p1y, p2y)) and
            (x <= max(p1x, p2x)) and (p1y!= p2y)):
            xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
            if p1x == p2x or x <= xints:
                inside = not inside
        p1x, p1y = p2x, p2y
    return inside
def is_edge_within_polygon(polygon, p1, p2, num_points=10):
    ##print(polygon[start], polygon[end])
    points = [(p1[0] + i*(p2[0]-p1[0])/num_points, p1[1] + i*(p2[1]-p1[1])/num_points) for i in range(1,num_points-1)]
    points_inside = [point_in_polygon(point, polygon) for point in points]
    return all(points_inside)
    #line_inside_polygon = check_line_in_polygon(line, polygon, num_points=10)

def distance(p1, p2, threshold=0.22):
    dist = math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
    return round(dist, 3)

def is_convex(polygon, idx):
    n = len(polygon)
    prev = polygon[idx - 1]
    curr = polygon[idx]
    next = polygon[(idx + 1) % n]
    cross_product = (curr[0] - prev[0]) * (next[1] - prev[1]) - (curr[1] - prev[1]) * (next[0] - prev[0])
    return cross_product > 0
def find_concave_nodes(polygon):
    return [i for i in range(len(polygon)) if not is_convex(polygon, i)]
def find_line_end_nodes(polygon):
    convex_nodes = [i for i in range(len(polygon)) if is_convex(polygon, i)]
    #print(convex_nodes)
    line_end_nodes = []

    if len(convex_nodes) < 2:
        return line_end_nodes

    current_set = []
    for i in range(len(convex_nodes) -1):
        current_node = convex_nodes[i]
        next_node = convex_nodes[i + 1]
        if next_node == current_node + 1 or (current_node == len(polygon) - 1 and next_node == 0):
            current_set.append(current_node)
            # #print("unluckily")
            # #print(next_node , current_node + 1, current_node, len(polygon) - 1 ,next_node )
        else:
            if len(current_set) >= 1:
                current_set.append(current_node)
                # Process the current set
                if len(current_set) > 1:
                    min_distance = float('inf')
                    line_end_node = None
                    for j in range(len(current_set) - 1):
                        dist = distance(polygon[current_set[j]], polygon[current_set[j + 1]])
                        #print("dist", dist)
                        #print("min_distance", min_distance)
                        if dist < min_distance:
                            min_distance = dist
                            line_end_node = current_set[j]
                            #print("line_end_node", line_end_node)
                    if line_end_node is not None:
                        line_end_nodes.append(line_end_node)
                        #print(line_end_nodes)
            current_set = []
        if current_node + 1 == len(polygon)-1:
            if len(current_set) >= 1:
                current_set.append(current_node)
                # Process the current set
                if len(current_set) > 1:
                    min_distance = float('inf')
                    line_end_node = None
                    for j in range(len(current_set) - 1):
                        dist = distance(polygon[current_set[j]], polygon[current_set[j + 1]])
                        if dist < min_distance:
                            min_distance = dist
                            line_end_node = current_set[j]
                    if line_end_node is not None:
                        line_end_nodes.append(line_end_node)
                        #print(line_end_nodes)
            current_set = []
    return line_end_nodes
def find_concave_nodes(polygon):
    return [i for i in range(len(polygon)) if not is_convex(polygon, i)]

class Graph:
    def __init__(self):
        self.graph = defaultdict(list)
        self.is_line_end = defaultdict(bool)

    def add_edge(self, u, v, length):
        self.graph[u].append((v, length))
        self.graph[v].append((u, length))

    def set_line_end(self, node, is_end):
        self.is_line_end[node] = is_end

def dijkstra(graph, start):
    distances = {node: float('inf') for node in graph}
    distances[start] = 0
    priority_queue = [(0, start)]
    shortest_paths = {node: [] for node in graph}
    shortest_paths[start] = [start]

    while priority_queue:
        current_distance, current_node = heapq.heappop(priority_queue)

        if current_distance > distances[current_node]:
            continue

        for neighbor, weight in graph[current_node]:
            distance = current_distance + weight

            if distance < distances[neighbor]:
                distances[neighbor] = distance
                shortest_paths[neighbor] = shortest_paths[current_node] + [neighbor]
                heapq.heappush(priority_queue, (distance, neighbor))

    return distances, shortest_paths

def longest_shortest_path(graph, is_line_end):
    line_end_nodes = [node for node, is_end in is_line_end.items() if is_end]
    max_length_info = [0, []]

    for i in range(len(line_end_nodes)):
        start = line_end_nodes[i]
        distances, shortest_paths = dijkstra(graph, start)

        for j in range(i + 1, len(line_end_nodes)):
            end = line_end_nodes[j]
            if distances[end] != float('inf') and distances[end] > max_length_info[0]:
                max_length_info[0] = round(distances[end],3)
                max_length_info[1] = shortest_paths[end]

    return max_length_info

def over_limit_distance_path(graph, is_line_end, limit_distance):
    line_end_nodes = [node for node, is_end in is_line_end.items() if is_end]
    result = []

    for i in range(len(line_end_nodes)):
        start = line_end_nodes[i]
        distances, shortest_paths = dijkstra(graph, start)

        for j in range(i + 1, len(line_end_nodes)):
            end = line_end_nodes[j]
            if distances[end]!= float('inf') and distances[end] > limit_distance:
                result.append([round(distances[end], 3), shortest_paths[end]])

    return result


def over_limit_distance_path(graph, is_line_end, limit_distance, Maxnumber):
    line_end_nodes = [node for node, is_end in is_line_end.items() if is_end]
    result = []

    for start in line_end_nodes:
        distances, shortest_paths = dijkstra(graph, start)

        for end in line_end_nodes:
            if start != end and distances[end] != float('inf') and distances[end] > limit_distance:
                heapq.heappush(result, (-distances[end], shortest_paths[end]))

    sorted_result = [heapq.heappop(result) for _ in range(len(result))]
    sorted_result.reverse()
    return sorted_result[:Maxnumber]


def generate_Metal_Vias(Metal_polygons, Via_Points):
    Metal_Vias = []
    visited = set()  
    for Metal_polygon in Metal_polygons:
        Metal_Via = []
        for index, Via_Point in enumerate(Via_Points):
            if Via_Point not in visited:
                if point_in_polygon(Via_Point, Metal_polygon):
                    Metal_Via.append(index)
                    visited.add(Via_Point)
        Metal_Vias.append(Metal_Via)
    return Metal_Vias

def generate_edges(Polygon, Polygon_prefix, Vias, Insidevias, Via_prefix):
    n = len(Polygon)
    m = len(Insidevias)
    edges = []
    for i in range(n):
        j = (i + 1) % n
        #edges.append((f"{Polygon_prefix}_{i}", f"{Polygon_prefix}_{j}", distance(Polygon[i], Polygon[j])))
        edges.append(("{0}_{1}".format(Polygon_prefix, i), "{0}_{1}".format(Polygon_prefix, j), distance(Polygon[i], Polygon[j])))
    for i in range(n):
        for j in range(i + 2, n):
            p1 = Polygon[i]
            p2 = Polygon[j]
            if is_edge_within_polygon(Polygon, p1, p2):
                #edges.append((f"{Polygon_prefix}_{i}", f"{Polygon_prefix}_{j}", round(dist, 3)(abs((p1[0]-(p2[0])))+abs((p1[1]-(p2[1]))))))
                edges.append(("{0}_{1}".format(Polygon_prefix, i), "{0}_{1}".format(Polygon_prefix, j), round(abs((p1[0] - (p2[0])))+abs((p1[1] - (p2[1]))),3)))

    for i in range(n):
        for j in range(m):
            p1 = Polygon[i]
            p2 = Vias[Insidevias[j]]
            if is_edge_within_polygon(Polygon, p1, p2):
                #edges.append((f"{Polygon_prefix}_{i}", f"{Via_prefix}_{Insidevias[j]}", distance(p1, p2)))
                edges.append(("{0}_{1}".format(Polygon_prefix, i), "{0}_{1}".format(Via_prefix, Insidevias[j]), distance(p1, p2)))
    return edges
def generate_edges(Polygon, Polygon_prefix, Vias, Insidevias, Via_prefix):
    edges = []
    #data = [x for point in Polygon for x in point]
    data = [x + random.uniform(-20, 20)/1011.01 for point in Polygon for x in point]
    triangles = earcut.earcut(data)
    for i in range(0, len(triangles), 3):
        a = triangles[i]
        b = triangles[i + 1]
        c = triangles[i + 2]
        edges.append(("{0}_{1}".format(Polygon_prefix, a), "{0}_{1}".format(Polygon_prefix, b), distance(Polygon[a], Polygon[b])))
        edges.append(("{0}_{1}".format(Polygon_prefix, a), "{0}_{1}".format(Polygon_prefix, c), distance(Polygon[a], Polygon[c])))
        edges.append(("{0}_{1}".format(Polygon_prefix, b), "{0}_{1}".format(Polygon_prefix, c), distance(Polygon[b], Polygon[c])))
    print(edges)
    return edges
def generate_connection_net(Metal_polygons, Metal_prefix, *args):
    Metal_edges = []
    Metal_via_list = []

    Via_Points_list = []
    Via_prefix_list = []
    
    for i in range(0, len(args), 2):
        Via_Points_list.append(args[i])
        Via_prefix_list.append(args[i + 1])
    for Via_Points in Via_Points_list:
        Metal_via_list.append(generate_Metal_Vias(Metal_polygons, Via_Points))
    
    for index, Metal_polygon in enumerate(Metal_polygons):
        print("start generate metal {}".format(index))
	print(len(Metal_polygon))
	print(Via_Points_list)
        for via_index in range(len(Via_Points_list)):
            Metal_edges.extend(
                generate_edges(
                    Metal_polygon, 
                    "{0}_{1}".format(Metal_prefix, index), 
                    Via_Points_list[via_index], 
                    Metal_via_list[via_index][index], 
                    "{0}_{1}".format(Via_prefix_list[via_index], 0)
                )
            )
    
    return Metal_edges

#find line_end for each group
def generate_line_end_nodes(Net_polygons, *All_metal_polygons_and_prefix):
    Net_line_ends = []
    visited = set()  
    for Net_polygon in Net_polygons:
        Net_line_end = []
        for metal_polygons, metal_prefix in All_metal_polygons_and_prefix:
            Metal_net_line_end = []
            for second_index, metal_polygon in enumerate(metal_polygons):
                for third_index, metal_point in enumerate(metal_polygon):
                    #if metal_point not in visited:
                    if tuple(metal_point) not in visited:
                        if is_convex(metal_polygon, third_index):
                            if point_in_polygon(metal_point, Net_polygon):
                                #Metal_net_line_end.append(f'{metal_prefix}_{second_index}_{third_index}')
                                Metal_net_line_end.append("{0}_{1}_{2}".format(metal_prefix, second_index, third_index))
                                visited.add(tuple(metal_point))
            Net_line_end.extend(Metal_net_line_end)
        Net_line_ends.append(Net_line_end)
    return Net_line_ends
def calculate_centers(polygons):
    centers = []
    for polygon in polygons:
        if len(polygon) == 4:  # Ensure it is a quadrilateral
            cx = sum(point[0] for point in polygon) / 4
            cy = sum(point[1] for point in polygon) / 4
            centers.append((cx, cy))
    return centers

def get_point_from_node(node, M1_polygons, M2_polygons, M3_polygons, M4_polygons, M5_polygons, V1_Points, V2_Points, V3_Points, V4_Points):
    parts = node.split('_')
    if parts[0] == "M1":
        poly_index = int(parts[1])
        point_index = int(parts[2])
        point = M1_polygons[poly_index][point_index]
    elif parts[0] == "M2":
        poly_index = int(parts[1])
        point_index = int(parts[2])
        point = M2_polygons[poly_index][point_index]
    elif parts[0] == "M3":
        poly_index = int(parts[1])
        point_index = int(parts[2])
        point = M3_polygons[poly_index][point_index]
    elif parts[0] == "M4":
        poly_index = int(parts[1])
        point_index = int(parts[2])
        point = M4_polygons[poly_index][point_index]
    elif parts[0] == "M5":
        poly_index = int(parts[1])
        point_index = int(parts[2])
        point = M5_polygons[poly_index][point_index]
    elif parts[0] == "V1":
        point_index = int(parts[2])
        point = V1_Points[point_index]
    elif parts[0] == "V2":
        point_index = int(parts[2])
        point = V2_Points[point_index]
    elif parts[0] == "V3":
        point_index = int(parts[2])
        point = V3_Points[point_index]
    elif parts[0] == "V4":
        point_index = int(parts[2])
        point = V4_Points[point_index]
    else:
        return None
    return (round(point[0], 3), round(point[1], 3))




#................................................................#
#Main 
Length_Limit = 2200
Length_Limit = 50
Area_Limit = 2200*0.08
Area_Limit = 500*0.08
Perimeter_Limit = 2200*2
Perimeter_Limit = 500*2
ground_names = [
    "VSS", "VSSCK", "VSSQ", "VG_BLP", "VSS_CLEAN", "SEL*", "AY*",
    "vss", "vssck", "vssq", "vg_blp", "VNR1", "VNR2", 
    "VNR3", "VNH", "VNRH1", "VNRH2", "VNSRH", "VPR1", 
    "VN", "VDD", "VDDQ", "VPSH", "VPRH", "VDSM1", 
    "VDSM2", "VDCK", "VDS?", "VPH", "VCC", "VDD_CLEAN", 
    "vdd", "vddq", "vpsh", "vprh", "vdsm1", "vdsm2", 
    "vdck", "vds?", "vph", "vcc", "VPSH_ASW1", "VPRH1", 
    "VPRH2", "A2C_VG_BLP", "VPR2"
]
tem_dir = r"/workspace/projs/dth/proj_lib/L13450163/check/run_dicectory/0722/"

#generate DRC_rule
checkbuffer_rule_name = r"{}buffer_check_debug".format(tem_dir)
generate_checkbuffer_rule(Length_Limit, Area_Limit, Perimeter_Limit, ground_names, checkbuffer_rule_name)
#generate DRC_cal
#layout_path = r"/workspace/projs/dth/proj_lib/L13450163/check/run_dicectory/0722/CORE_TILE_PERI.gds"
layout_path = r"/workspace/projs/dth/proj_lib/L13450163/check/run_dicectory/0722/DTH_PERI_RevC.gds"
#layout_path = r"/workspace/projs/dth/proj_lib/L13450163/check/run_dicectory/0722/temp.gds"
layout_primary = r"DTH_PERI_RevC"
#layout_primary = r"CORE_TILE_PERI"
#layout_primary = r"temp"
drc_results_database = r"DTH_PERI_RevC.drc.results"
#drc_results_database = r"CORE_TILE_PERI.drc.results"
#drc_results_database = r"temp.drc.results"
cal_file_name = r"{}drc_cal_".format(tem_dir)
output_drc_filename = r"{}_{}buffer_skill.drc.results".format(tem_dir,layout_primary)
generate_drc_cal(layout_path, layout_primary, drc_results_database, checkbuffer_rule_name, cal_file_name)
os.chdir(tem_dir)
#os.system("calibre -drc -hier -hyper -turbo -turbo_all {}".format(cal_file_name))
print("start read data!")
data = glue.parse_document(drc_results_database)
print("read data finished!")
variables = glue.output_polygons(data)
#M1_polygons = variables.get('M1_polygons')
#M2_polygons = variables.get('M2_polygons')
#M3_polygons = variables.get('M2_polygons')
#M4_polygons = variables.get('M2_polygons')
#M5_polygons = variables.get('M2_polygons')
#Net_polygons = variables.get('Net_polygons')
#V1_polygons = variables.get('V1C_polygons') + variables.get('V1P_polygons')
#V1_Points = calculate_centers(V1_polygons)
#V2_polygons = variables.get('V2C_polygons') + variables.get('V2P_polygons')
#V2_Points = calculate_centers(V2_polygons)
#V3_polygons = variables.get('V3C_polygons') + variables.get('V3P_polygons')
#V3_Points = calculate_centers(V3_polygons)
#V4_polygons = variables.get('V4C_polygons') + variables.get('V4P_polygons')
#V4_Points = calculate_centers(V4_polygons)
print("create variables finished!")
#print(V1_Points )	
Net_polygons = variables.get('Net_polygons')
Net_polygons = [polygon for polygon in Net_polygons if len(polygon) < 5000]
Net_polygons = simplify.simplify_polygons(Net_polygons,0.06)
Net_vias = []
print("start skellize")
Net_polygons_skels = []
for Net_polygon in Net_polygons:
    print(Net_polygon)
    Net_polygons_skel = polyskel.skeletonize(Net_polygon, Net_vias)
    Net_polygons_skels.append(Net_polygons_skel)
print("skellize_finished")
#print(len(Net_polygons[0]))
M1_polygons = Net_polygons
#M1_polygons = [polygon for polygon in Net_polygons if len(polygon) < 5000]
M2_polygons = []
M3_polygons = []
M4_polygons = []
M5_polygons = []
V1_polygons = []
V1_Points = calculate_centers(V1_polygons)
V2_polygons = []
V2_Points = calculate_centers(V2_polygons)
V3_polygons = []
V3_Points = calculate_centers(V3_polygons)
V4_polygons = []
V4_Points = calculate_centers(V4_polygons)
output_polygon_filename = r"/home/you_pan/Desktop/polygon.txt"
#with open(output_polygon_filename, 'w') as file:
    #file.write(str(Net_polygons[0]))
    #file.write(str(Net_polygons[1]))
    #file.write(str(Net_polygons[2]))
    #file.write(str(Net_polygons[3]))
print("export polygons finished")
M1_connection_net = generate_connection_net(M1_polygons, "M1", V1_Points, "V1")
M2_connection_net = generate_connection_net(M2_polygons, "M2", V1_Points, "V1", V2_Points, "V2")
M3_connection_net = generate_connection_net(M3_polygons, "M3", V2_Points, "V2", V3_Points, "V3")
M4_connection_net = generate_connection_net(M4_polygons, "M4", V3_Points, "V3", V4_Points, "V4")
M5_connection_net = generate_connection_net(M5_polygons, "M5", V4_Points, "V4")
total_connection_net = M1_connection_net + M2_connection_net + M3_connection_net + M4_connection_net + M5_connection_net
print("generate connections of metals nodes finished!")
total_line_end_nodes = generate_line_end_nodes(Net_polygons, [M1_polygons,"M1"], [M2_polygons,"M2"], [M3_polygons,"M3"], [M4_polygons,"M4"], [M5_polygons,"M5"])
print("make connections between nodes finished!")
g = Graph()
#for u, v, length in total_connection_net:
    #g.add_edge(u, v, length)



#with open(drc_filename, 'w') as file:
#     file.write("{} 10000\n".format("buffer_skill"))
#     for net_index, total_line_end_node in enumerate(total_line_end_nodes):
#         g = Graph()
#         for u, v, length in total_connection_net:
#             g.add_edge(u, v, length)
#         for node in total_line_end_node:
#             g.set_line_end(node, True)
#         max_length_info = longest_shortest_path(g.graph, g.is_line_end)
#         over_limit_distance_info = over_limit_distance_path(g.graph, g.is_line_end, Length_Limit)
#         longest_path_nodes = max_length_info[1] 
#         longest_path_points = [get_point_from_node(node, M1_polygons, M2_polygons, M3_polygons, M4_polygons, M5_polygons, V1_Points, V2_Points, V3_Points, V4_Points) for node in longest_path_nodes]
# 	#print(longest_path_points)
#         longest_path_polygons = export_to_drc.generate_path_points(longest_path_points, 0.02, decimals=3)
#         #print("Longest path length:", max_length_info[0])
#         #print("Longest path nodes:", max_length_info[1])
#         #print("Longest path points:", longest_path_points)
#         #print("Longest path points:", longest_path_polygons)
#         file.write("Max_Length_{}_net_{}\n".format(max_length_info[0], net_index))
#         file.write("1 1 1 Jun 28 09:25:36 2024 \n")
#         if longest_path_polygons is not None:
#                 file.write("'net_{}' {{{}}}\n".format(net_index,"Max Length = {}".format(max_length_info[0])))
#                 file.write("p 1 {}\n".format(len(longest_path_polygons)))
#                 for point in longest_path_polygons:
#                     file.write("{} {}\n".format(int(point[0]*10000),int(point[1]*10000)))

with open(output_drc_filename, 'w') as file:
    file.write("{} 10000\n".format(layout_primary))
    for net_index, total_line_end_node in enumerate(total_line_end_nodes):
        g = Graph()
        for u, v, length in total_connection_net:
            g.add_edge(u, v, length)
        for node in total_line_end_node:
            g.set_line_end(node, True)
        max_length_info = longest_shortest_path(g.graph, g.is_line_end)
        over_limit_distance_info = over_limit_distance_path(g.graph, g.is_line_end, Length_Limit,10)
        over_limit_distance_info = over_limit_distance_path(g.graph, g.is_line_end, 1,100)
        for over_distance, over_path_nodes in over_limit_distance_info:
            #longest_path_nodes = max_length_info[1] 
            longest_path_points = [get_point_from_node(node, M1_polygons, M2_polygons, M3_polygons, M4_polygons, M5_polygons, V1_Points, V2_Points, V3_Points, V4_Points) for node in over_path_nodes]
        #print(longest_path_points)
            longest_path_polygons = export_to_drc.generate_path_points(longest_path_points, 0.02, decimals=3)
            #print("Longest path length:", max_length_info[0])
            #print("Longest path nodes:", max_length_info[1])
            #print("Longest path points:", longest_path_points)
            #print("Longest path points:", longest_path_polygons)
            file.write("Max_Length_{}_net_{}\n".format(over_distance, net_index))
            file.write("1 1 1 Jun 28 09:25:36 2024 \n")
            if longest_path_polygons is not None:
                    file.write("'net_{}' {{{}}}\n".format(net_index,"Max Length = {}".format(over_distance)))
                    file.write("p 1 {}\n".format(len(longest_path_polygons)))
                    for point in longest_path_polygons:
                        file.write("{} {}\n".format(int(point[0]*10000),int(point[1]*10000)))
