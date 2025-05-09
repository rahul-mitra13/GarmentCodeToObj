import argparse
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
import openmesh as om
import triangle
import polyscope as ps

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Read a JSON file and maximum triangle area from the command line.")
    parser.add_argument('file', help="Path to the JSON file")
    parser.add_argument('a', help="Maximum traingle area")

    # Parse the arguments
    args = parser.parse_args()

    #store the user-specified max area 
    max_area = args.a

    # Read and parse the JSON file
    try:
        with open(args.file, 'r') as f:
            data = json.load(f)
            panels = data["pattern"]["panels"]
            stitches = data["pattern"]["stitches"]
            #make a mesh from the panels
            makeMeshes(panels, stitches, max_area)
        
        # Print the parsed JSON data
        print("Successfully read JSON file:")

    except FileNotFoundError:
        print(f"Error: File '{args.file}' not found.")
    except json.JSONDecodeError:
        print(f"Error: '{args.file}' is not a valid JSON file.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

#make a number of meshes from panels so that they agree across stitched edges
def makeMeshes(panels, stitches, max_area):
    #create a new dictionary where we add panel vertices from edges (this is the correct ordering we want)
    #the idea is to find the correct ordering of vertices and re-order vertices when you make them from stitches
    panel_vertices_from_edges = {}
    #segments we want to retain in the triangulation
    panel_segments_from_edges = {}
    #create a new dictionary where we will map vertices to panels (these panel vertices come from stitches)
    panel_vertices_from_stitches = {}
    #create a new dicitionary where we will store vertex mappings between panels
    vertex_mappings = {}
    #create a new dictionary where we will store panels to triangulations (this is needed to create the global mesh)
    #these are triangulations from the triangle library
    panel_triangulations = {}    
    #min and max number of stitches across an edge 
    min_stitches_on_edge = 10
    max_stitches_on_edge = 100
    #min and max edge lengths in the original garment 
    min_edge_length = 10000
    max_edge_length = 0
    #first set panel_vertices_from_edges to an empty list 
    #also set panel_segments_from_edges to an empty list
    for stitch in stitches:
        panel1 = stitch[0]['panel']
        panel2 = stitch[1]['panel']
        panel_vertices_from_edges.setdefault(panel1, [])
        panel_vertices_from_edges.setdefault(panel2, [])
        panel_segments_from_edges.setdefault(panel1, [])
        panel_segments_from_edges.setdefault(panel2, [])

    #first find min and max edge lengths for normalization 
    for panel in panels:
        vertices = panels[panel]["vertices"]
        edges = panels[panel]["edges"]
        vertices_to_triangulate = []
        segments = []
        for e in edges: 
            #handle the case when the edge is a quadratic bezier curve
            if "curvature" in e:
                v1 = vertices[e["endpoints"][0]]
                v2 = vertices[e["endpoints"][1]]
                #the curvature coords
                #treating the edge as a vector (1, 0) according to the garment code specification
                curvature_coords = e["curvature"]
                #v3 is our control point
                v3 = [None] * 2
                v3 = relative_to_abs_coord(v1, v2, curvature_coords)
                length = bezier_arc_length(v1, v3, v2)
                if (length < min_edge_length):
                    min_edge_length = length
                if (length > max_edge_length):
                    max_edge_length = length    
            else:
                v1 = vertices[e["endpoints"][0]]
                v2 = vertices[e["endpoints"][1]]
                length = line_length(v1, v2)
                if (length < min_edge_length):
                    min_edge_length = length
                if (length > max_edge_length):
                    max_edge_length = length 
    
    
    #now make vertices based on edge ordering
    for panel in panels:
        vertices = panels[panel]["vertices"]
        edges = panels[panel]["edges"]
        vertices_to_triangulate = []
        segments = []
        for e in edges: 
            #handle the case when the edge is a quadratic bezier curve
            if "curvature" in e:
                v1 = vertices[e["endpoints"][0]]
                v2 = vertices[e["endpoints"][1]]
                #the curvature coords
                #treating the edge as a vector (1, 0) according to the garment code specification
                curvature_coords = e["curvature"]
                #v3 is our control point
                v3 = [None] * 2
                v3 = relative_to_abs_coord(v1, v2, curvature_coords)
                length = bezier_arc_length(v1, v3, v2)
                #print("num vertices on stitch: ", map_to_int_range(length, min_edge_length, max_edge_length, min_stitches_on_edge, max_stitches_on_edge))
                num_vertices_on_stitch = map_to_int_range(length, min_edge_length, max_edge_length, min_stitches_on_edge, max_stitches_on_edge)
                uniform_distance_points = sample_quadratic_bezier(v1, v3, v2, num_vertices_on_stitch)
                #push back the computed points on the Bezier curve to the vertices we want to triangulate 
                for p in uniform_distance_points:
                    panel_vertices_from_edges[panel].append(p)
            else:
                v1 = vertices[e["endpoints"][0]]
                v2 = vertices[e["endpoints"][1]]
                length = line_length(v1, v2)
                #print("num vertices on stitch: ", map_to_int_range(length, min_edge_length, max_edge_length, min_stitches_on_edge, max_stitches_on_edge))
                num_vertices_on_stitch = map_to_int_range(length, min_edge_length, max_edge_length, min_stitches_on_edge, max_stitches_on_edge)
                uniform_distance_points = sample_equally_spaced_points_on_line(v1, v2, num_vertices_on_stitch)
                for p in uniform_distance_points:
                    panel_vertices_from_edges[panel].append(p)                
        #remove co-located vertices (probably not needed?)
        panel_vertices_from_edges[panel] = remove_colocated_vertices(panel_vertices_from_edges[panel])
        #make the segments for a specific panel 
        for i in range(len(panel_vertices_from_edges[panel])):
            panel_segments_from_edges[panel].append((i, (i + 1) % len(panel_vertices_from_edges[panel])))
    
    #first initialize all the panel vertices to an empty list 
    #this is a bit redundant but okay for now
    for stitch in stitches:
        panel1 = stitch[0]['panel']
        panel2 = stitch[1]['panel']
        panel_vertices_from_stitches.setdefault(panel1, [])
        panel_vertices_from_stitches.setdefault(panel2, [])
        vertex_mappings.setdefault(panel1, set())
        vertex_mappings.setdefault(panel2, set())

    for stitch in stitches: 
        #the first edge and first panel in the stitch 
        panel1 = stitch[0]['panel']
        edge1 = panels[panel1]["edges"][stitch[0]["edge"]]
        #the second edge and second panel in the stitch 
        panel2 = stitch[1]['panel']
        edge2 = panels[panel2]["edges"][stitch[1]["edge"]]
        #no matter what the stitch is, map the endpoints first
        #first panel
        panel1vertex1 = panels[panel1]["vertices"][edge1["endpoints"][0]]
        panel1vertex2 = panels[panel1]["vertices"][edge1["endpoints"][1]]
        #second panel
        panel2vertex1 = panels[panel2]["vertices"][edge2["endpoints"][0]]
        panel2vertex2 = panels[panel2]["vertices"][edge2["endpoints"][1]]
        
        index_1_from_panel_1 = find_colocated_vertex_index(panel1vertex1, panel_vertices_from_edges[panel1])
        index_2_from_panel_1 = find_colocated_vertex_index(panel1vertex2, panel_vertices_from_edges[panel1])
        index_1_from_panel_2 = find_colocated_vertex_index(panel2vertex1, panel_vertices_from_edges[panel2])
        index_2_from_panel_2 = find_colocated_vertex_index(panel2vertex2, panel_vertices_from_edges[panel2])

        vertex_mappings[panel1].add((panel2, (index_1_from_panel_1, index_2_from_panel_2)))
        vertex_mappings[panel1].add((panel2, (index_2_from_panel_1, index_1_from_panel_2)))
        vertex_mappings[panel2].add((panel1, (index_1_from_panel_2, index_2_from_panel_1)))
        vertex_mappings[panel2].add((panel1, (index_2_from_panel_2, index_1_from_panel_1)))
        
        
        #both edges on the stitch are specified by quadratic bezier curves 
        if "curvature" in edge1 and "curvature" in edge2: 
            panel1CurvatureCoords = edge1["curvature"]
            #control point
            panel1ControlPoint = [None] * 2
            #this is treating the curvature coords as (1, 0) [GarmentCode specification]
            panel1ControlPoint = relative_to_abs_coord(panel1vertex1, panel1vertex2, panel1CurvatureCoords)
            #second panel
            panel2CurvatureCoords = edge2["curvature"]
            #control point
            panel2ControlPoint = [None] * 2
            #this is treating the curvature coords as (1, 0) [GarmentCode specification]
            panel2ControlPoint = relative_to_abs_coord(panel2vertex1, panel2vertex2, panel2CurvatureCoords)
            
            #just find the length of one of them since we want the length to be equal anyway
            length = bezier_arc_length(panel1vertex1, panel1ControlPoint, panel1vertex2)
            num_vertices_on_stitch = map_to_int_range(length, min_edge_length, max_edge_length, min_stitches_on_edge, max_stitches_on_edge)
            
            #generate vertices on the bezier curve
            #hopefully the triangle library maintains an ordering on vertices 
            uniform_distance_points_panel1 = sample_quadratic_bezier(panel1vertex1, panel1ControlPoint, panel1vertex2, num_vertices_on_stitch)
            #check if you need to swap to keep indexing correct
            #swap panel 2 vertices since they're oppositely ordered across a stitch
            panel2vertex1, panel2vertex2 = panel2vertex2, panel2vertex1
            uniform_distance_points_panel2 = sample_quadratic_bezier(panel2vertex1, panel2ControlPoint, panel2vertex2, num_vertices_on_stitch)
            #push back the computed points on the Bezier and find mapping 
            for i in range(len(uniform_distance_points_panel1)):
                panel_vertices_from_stitches[panel1].append(uniform_distance_points_panel1[i])
                index_from_panel_1 = find_colocated_vertex_index(uniform_distance_points_panel1[i], panel_vertices_from_edges[panel1])
                panel_vertices_from_stitches[panel2].append(uniform_distance_points_panel2[i])
                index_from_panel_2 = find_colocated_vertex_index(uniform_distance_points_panel2[i], panel_vertices_from_edges[panel2])
                vertex_mappings[panel1].add((panel2, (index_from_panel_1, index_from_panel_2)))
                vertex_mappings[panel2].add((panel1, (index_from_panel_2, index_from_panel_1)))

        #both edges on the stitch are specified by a straight curve 
        elif "curvature" not in edge1 and "curvature" not in edge2: 

            length = line_length(panel1vertex1, panel1vertex2)
            num_vertices_on_stitch = map_to_int_range(length, min_edge_length, max_edge_length, min_stitches_on_edge, max_stitches_on_edge)

            #generate vertices on the straight edge 
            uniform_distance_points_panel1 = sample_equally_spaced_points_on_line(panel1vertex1, panel1vertex2, num_vertices_on_stitch)
            #swap panel 2 vertices since they're oppositely ordered across a stitch
            panel2vertex1, panel2vertex2 = panel2vertex2, panel2vertex1
            uniform_distance_points_panel2 = sample_equally_spaced_points_on_line(panel2vertex1, panel2vertex2, num_vertices_on_stitch)
            #push back the computed points on the straight line and find mapping
            for i in range(len(uniform_distance_points_panel1)):
                panel_vertices_from_stitches[panel1].append(uniform_distance_points_panel1[i])
                index_from_panel_1 = find_colocated_vertex_index(uniform_distance_points_panel1[i], panel_vertices_from_edges[panel1])
                panel_vertices_from_stitches[panel2].append(uniform_distance_points_panel2[i])
                index_from_panel_2 = find_colocated_vertex_index(uniform_distance_points_panel2[i], panel_vertices_from_edges[panel2])
                vertex_mappings[panel1].add((panel2, (index_from_panel_1, index_from_panel_2)))
                vertex_mappings[panel2].add((panel1, (index_from_panel_2, index_from_panel_1)))

        #edge 1 on the stitch is specified by a straight curve and edge 2 is a quadratic bezier curve
        elif "curvature" not in edge1 and "curvature" in edge2: 
            #don't think this can ever happen but just use length of straight line? 
            length = line_length(panel1vertex1, panel1vertex2)
            num_vertices_on_stitch = map_to_int_range(length, min_edge_length, max_edge_length, min_stitches_on_edge, max_stitches_on_edge)
            uniform_distance_points_panel1 = sample_equally_spaced_points_on_line(panel1vertex1, panel1vertex2, num_vertices_on_stitch)
            panel2ControlPoint = [None] * 2
            #this is treating the curvature coords as (1, 0) [GarmentCode specification]
            panel2ControlPoint = relative_to_abs_coord(panel2vertex1, panel2vertex2, panel2CurvatureCoords)
            #swap panel 2 vertices since they're oppositely ordered across a stitch
            panel2vertex1, panel2vertex2 = panel2vertex2, panel2vertex1
            uniform_distance_points_panel2 = sample_quadratic_bezier(panel2vertex1, panel2ControlPoint, panel2vertex2, num_vertices_on_stitch)
            #push back the computed points and find the mappings
            for i in range(len(uniform_distance_points_panel1)):
                panel_vertices_from_stitches[panel1].append(uniform_distance_points_panel1[i])
                index_from_panel_1 = find_colocated_vertex_index(uniform_distance_points_panel1[i], panel_vertices_from_edges[panel1])
                panel_vertices_from_stitches[panel2].append(uniform_distance_points_panel2[i])
                index_from_panel_2 = find_colocated_vertex_index(uniform_distance_points_panel2[i], panel_vertices_from_edges[panel2])
                vertex_mappings[panel1].add((panel2, (index_from_panel_1, index_from_panel_2)))
                vertex_mappings[panel2].add((panel1, (index_from_panel_2, index_from_panel_1)))
        
        #edge 1 on the stitch is specified by a quadratic bezier curve and edge 2 is specified by a straight curve
        #I don't think this can ever happen 
        elif "curvature" in edge1 and "curvature" not in edge2:
            #don't think this can ever happen but just use length of straight line? 
            length = line_length(panel2vertex1, panel2vertex2)
            num_vertices_on_stitch = map_to_int_range(length, min_edge_length, max_edge_length, min_stitches_on_edge, max_stitches_on_edge)
            panel1CurvatureCoords = edge1["curvature"]
            #control point
            panel1ControlPoint = [None] * 2
            #this is treating the curvature coords as (1, 0) [GarmentCode specification]
            panel1ControlPoint = relative_to_abs_coord(panel1vertex1, panel1vertex2, panel1CurvatureCoords)
            uniform_distance_points_panel1 = sample_quadratic_bezier(panel1vertex1, panel1ControlPoint, panel1vertex2, num_vertices_on_stitch)
            #swap panel 2 vertices since they're oppositely ordered across a stitch
            panel2vertex1, panel2vertex2 = panel2vertex2, panel2vertex1
            uniform_distance_points_panel2 = sample_equally_spaced_points_on_line(panel2vertex1, panel2vertex2, num_vertices_on_stitch)
            #push back the computed points and find the mappings
            for i in range(len(uniform_distance_points_panel1)):
                panel_vertices_from_stitches[panel1].append(uniform_distance_points_panel1[i])
                index_from_panel_1 = find_colocated_vertex_index(uniform_distance_points_panel1[i], panel_vertices_from_edges[panel1])
                panel_vertices_from_stitches[panel2].append(uniform_distance_points_panel2[i])
                index_from_panel_2 = find_colocated_vertex_index(uniform_distance_points_panel2[i], panel_vertices_from_edges[panel2])
                vertex_mappings[panel1].add((panel2, (index_from_panel_1, index_from_panel_2)))
                vertex_mappings[panel2].add((panel1, (index_from_panel_2, index_from_panel_1)))
                
    #triangulate the panels using CDT
    for key in panel_vertices_from_edges:
        A = dict(vertices=panel_vertices_from_edges[key], segments=panel_segments_from_edges[key])
        print("Triangulating " + key + "....")
        # Perform constrained Delaunay triangulation
        #can try and improve the quality of the triangulation here
        B = triangle.triangulate(A, 'peq1.0a' + str(max_area) + 'C')
        panel_triangulations[key] = B
        #read the translation and rotation to construct global mesh
        translation = np.array(panels[key]["translation"])
        rotation = np.array(panels[key]["rotation"])
        #debug
        # print("Running some sanity checks...")
        # print("Does the triangulation have duplicate vertices? ", has_duplicate_vertices(B['vertices']))
        # print("Does the triangulation have duplicate triangles? ", has_duplicate_faces(B['triangles']))
        # print("Does the triangulation have duplicate edges? ", has_duplicate_edges(B['edges']))
    
    #make the global obj
    makeGlobalObj(panel_triangulations, vertex_mappings, panels)

#-----------HELPER FUNCTION DEFINITIONS-------------------#
#remove co-located vertices to clean up input 
def remove_colocated_vertices(points, tolerance=1e-8):
    # Convert to numpy array if not already
    points = np.array(points)
    
    # Use np.unique with a small tolerance
    _, unique_indices = np.unique(points.round(decimals=int(-np.log10(tolerance))), axis=0, return_index=True)
    
    # Return the unique points
    return points[np.sort(unique_indices)]


def sample_equally_spaced_points_on_line(start_point, end_point, num_points):
    """
    Sample equally spaced points on a line between start_point and end_point. This includes first and last point

    Parameters:
    - start_point: Tuple (x1, y1) representing the start point of the line.
    - end_point: Tuple (x2, y2) representing the end point of the line.
    - num_points: Number of equally spaced points to sample, including the endpoints.

    Returns:
    - points: A list of tuples representing the sampled points.
    """
    x1, y1 = start_point
    x2, y2 = end_point
    
    # Generate equally spaced points in x and y
    x_points = np.linspace(x1, x2, num_points)
    y_points = np.linspace(y1, y2, num_points)
    
    # Combine x and y points into a list of tuples
    points = list(zip(x_points, y_points))
    
    return points

def find_colocated_vertex_index(vertex, vertex_list, tol=1e-8):
    """
    Find the index of a co-located vertex in a list of vertices considering floating-point errors.

    Parameters:
    - vertex: A tuple representing the vertex to search for (e.g., (x, y)).
    - vertex_list: A list of tuples representing the list of vertices.
    - tol: Tolerance level for floating-point comparison (default is 1e-9).

    Returns:
    - index: The index of the co-located vertex if found, otherwise None.
    """
    for i, v in enumerate(vertex_list):
        if all(abs(a - b) < tol for a, b in zip(v, vertex)):
            return i
    return None

def are_vertices_colocated(vertex1, vertex2, tol=1e-9):
    """
    Check whether two vertices are co-located in space within a given tolerance.

    Parameters:
    - vertex1: A tuple (x, y) representing the first vertex.
    - vertex2: A tuple (x, y) representing the second vertex.
    - tol: A float representing the tolerance for floating-point comparison.

    Returns:
    - bool: True if the vertices are co-located within the tolerance, False otherwise.
    """
    vertex1 = np.array(vertex1)
    vertex2 = np.array(vertex2)
    return np.linalg.norm(vertex1 - vertex2) < tol

def relative_to_abs_coord(start, end, curvature_coords):
    """
    Derives absolute coordinates of Bezier control point given curvature coordinates in local edge frame
    See specifaction of GarmentCode (https://github.com/maria-korosteleva/GarmentCode?tab=readme-ov-file) for more details
    """
    start = np.array(start)
    end = np.array(end)
    curvature_coords = np.array(curvature_coords)
    edge = end - start
    edge_perp = np.array([-edge[1], edge[0]])

    control_start = start + curvature_coords[0] * edge
    control_point = control_start + curvature_coords[1] * edge_perp

    return control_point

#-------------QUADRATIC BEZIER STUFF---------------------#
def quadratic_bezier(t, p0, p1, p2):
    """
    Compute the coordinates of a point on a quadratic Bezier curve at parameter t.

    Parameters:
    - t: The parameter (0 <= t <= 1)
    - p0, p1, p2: Control points of the quadratic Bezier curve as tuples (x, y)

    Returns:
    - (x, y): The coordinates of the point on the Bezier curve at parameter t.
    """
    x = (1 - t)**2 * p0[0] + 2 * (1 - t) * t * p1[0] + t**2 * p2[0]
    y = (1 - t)**2 * p0[1] + 2 * (1 - t) * t * p1[1] + t**2 * p2[1]
    return x, y

def bezier_arc_length(p0, p1, p2):
    """
    Compute the arc length of the quadratic Bezier curve.

    Parameters:
    - p0, p1, p2: Control points of the quadratic Bezier curve as tuples (x, y)

    Returns:
    - length: The arc length of the Bezier curve.
    """
    def integrand(t):
        dx_dt = 2 * (1 - t) * (p1[0] - p0[0]) + 2 * t * (p2[0] - p1[0])
        dy_dt = 2 * (1 - t) * (p1[1] - p0[1]) + 2 * t * (p2[1] - p1[1])
        return np.sqrt(dx_dt**2 + dy_dt**2)
    
    length, _ = quad(integrand, 0, 1)
    return length

def sample_quadratic_bezier(p0, p1, p2, num_points):
    """
    Sample n equally spaced points on a quadratic Bezier curve. This includes first and last point

    Parameters:
    - p0, p1, p2: Control points of the quadratic Bezier curve as tuples (x, y)
    - num_points: Number of equally spaced points to sample

    Returns:
    - samples: List of tuples representing the sampled points [(x, y), (x, y), ...]
    """
    length = bezier_arc_length(p0, p1, p2)
    segment_length = length / (num_points - 1)

    def arc_length_parameter(t, target_length):
        def integrand(t):
            dx_dt = 2 * (1 - t) * (p1[0] - p0[0]) + 2 * t * (p2[0] - p1[0])
            dy_dt = 2 * (1 - t) * (p1[1] - p0[1]) + 2 * t * (p2[1] - p1[1])
            return np.sqrt(dx_dt**2 + dy_dt**2)
        current_length, _ = quad(integrand, 0, t)
        return current_length - target_length

    samples = [p0]
    for i in range(1, num_points - 1):
        target_length = i * segment_length
        t = brentq(arc_length_parameter, 0, 1, args=(target_length,))
        samples.append(quadratic_bezier(t, p0, p1, p2))
    samples.append(p2)

    return samples


def line_length(v1, v2):
    """
    Find the length of a straight line
    """
    point1 = np.array(v1)
    point2 = np.array(v2)
    length = np.linalg.norm(point2 - point1)
    return length

def map_to_int_range(value, original_min, original_max, new_min, new_max):
    """
    Takes in a value and maps it to a new range (integer)
    """
    # Normalize the value to the range [0, 1]
    normalized_value = (value - original_min) / (original_max - original_min)
    
    # Scale and shift the value to be within the new integer range [new_min, new_max]
    scaled_value = normalized_value * (new_max - new_min) + new_min
    
    # Round the result to the nearest integer
    integer_value = round(scaled_value)
    
    return integer_value

def make_halfedge_datastructure(vertices, faces):
    """
    Takes in a list of vertices and faces and creates a halfedge mesh 
    using openmesh
    """
    # Create an empty TriMesh
    mesh = om.TriMesh()
    # Add vertices to the mesh
    vertex_handles = []
    for vertex in vertices:
        vh = mesh.add_vertex(vertex)
        vertex_handles.append(vh)
 
    # Add faces to the mesh using the vertex handles
    for face in faces:
        mesh.add_face([vertex_handles[vid] for vid in face])
    
    print("Vertices:", mesh.n_vertices())
    print("Edges:", mesh.n_edges())
    print("Faces:", mesh.n_faces())

def makeGlobalObj(panel_triangulations, vertex_mappings, panels):
    """
    This function makes a global obj. Subjects the vertices to the global rotation and translation specified in the 
    GarmentCode JSON file 
    """

    #first map local vertex indices to global indices 
    global_index = 0
    local_to_global_index_map = {}
    for key in panel_triangulations:
        for i in range(len(panel_triangulations[key]['vertices'])):
            #this doesn't work because one vertex from a single panel might get mapped to multiple panels
            #^actually I think the above is fine (have to convince myself)
            #such as the middle vertex from a torso to two pants panel
            local_to_global_index_map[(key, i)] = global_index
            global_index += 1

    file = open('global.obj', 'w')
    vertices = []
    faces = []
    num_total_vertices = 0
    for key in panel_triangulations:
        #read the translation and rotation to construct global mesh
        translation = np.array(panels[key]["translation"])
        rotation = np.array(panels[key]["rotation"])
        # Write vertices
        for point in panel_triangulations[key]['vertices']:
            r = R.from_euler('xyz', rotation, degrees=True)
            updated_point = np.array([point[0], point[1], 0])
            updated_point = r.apply(updated_point)
            updated_point = updated_point + translation
            vertices.append(updated_point)
            file.write(f"v {updated_point[0]} {updated_point[1]} {updated_point[2]}\n")
        # Write faces
        for simplex in panel_triangulations[key]['triangles']:
            faces.append([simplex[0]+num_total_vertices, simplex[1]+num_total_vertices, simplex[2]+num_total_vertices])
            file.write(f"f {simplex[0]+1+num_total_vertices} {simplex[1]+1+num_total_vertices} {simplex[2]+1+num_total_vertices}\n")
        num_total_vertices += len(panel_triangulations[key]['vertices'])
    
    vertices = np.asarray(vertices)

    global_vertex_mappings = []
    #update the vertex mappings to contain global indices (don't create bi-directional mappings)
    used_keys = set()
    for key in vertex_mappings:
        for mapping in vertex_mappings[key]:
            #if something has showed up as a key ensure it doesn't show up as a value (this prevents something like (1 -> 15) and (15 -> 1) showing up)
            #need to add something so that we can have multiple duplicate keys 
            if (local_to_global_index_map[(mapping[0], mapping[1][1])]) not in used_keys:
                global_vertex_mappings.append((local_to_global_index_map[(key, mapping[1][0])], local_to_global_index_map[(mapping[0], mapping[1][1])]))
                used_keys.add(local_to_global_index_map[(key, mapping[1][0])])
    
    #debug
    # print("global vertex mapping:")
    # print(global_vertex_mappings)
    # print("Number of global vertex mappings: ", len(global_vertex_mappings))
    # print("Total number of vertices: ", num_total_vertices)

    # with open('global_vertex_mappings.txt', 'w') as f:
    #     for tup in global_vertex_mappings:
    #         f.write(f"{tup}\n")

    for tup in global_vertex_mappings:
        file.write(f"m {tup}\n")
    

    show_mesh(vertices, faces, global_vertex_mappings)

def has_duplicate_vertices(vertices, tol=1e-9):
    """
    Check for duplicate vertex positions in a list of vertex positions.
    
    Parameters:
    vertices (list of tuple): List of vertex positions, where each position is a tuple (x, y, z).
    tol (float): Tolerance for floating point comparison.
    
    Returns:
    bool: True if there are duplicate vertices, False otherwise.
    """
    seen = set()
    for vertex in vertices:
        # Round the vertex coordinates to handle floating-point precision issues
        rounded_vertex = tuple(np.round(vertex, decimals=int(-np.log10(tol))))
        if rounded_vertex in seen:
            return True
        seen.add(rounded_vertex)
    return False

def has_duplicate_faces(faces):
    """
    Check for duplicate faces in a list of faces.
    
    Parameters:
    faces (list of list or tuple): List of faces, where each face is a list or tuple of vertex indices.
    
    Returns:
    bool: True if there are duplicate faces, False otherwise.
    """
    seen = set()
    for face in faces:
        # Sort the vertex indices to handle different permutations of the same face
        sorted_face = tuple(sorted(face))
        if sorted_face in seen:
            return True
        seen.add(sorted_face)
    return False

def has_duplicate_edges(edges):
    """
    Check for duplicate edges in a list of edges.
    
    Parameters:
    edges (list of tuple): List of edges, where each edge is a tuple of vertex indices.
    
    Returns:
    bool: True if there are duplicate edges, False otherwise.
    """
    seen = set()
    for edge in edges:
        # Sort the vertex indices to handle undirected edges
        sorted_edge = tuple(sorted(edge))
        if sorted_edge in seen:
            return True
        seen.add(sorted_edge)
    return False

def show_mesh(vertices, faces, global_vertex_mappings):
    """
    Renders a mesh using polyscope and also renders the "virtual connections" 
    """
    ps.init()
    #register the surface mesh
    ps_mesh = ps.register_surface_mesh("global mesh", vertices, faces)
    #render "virtual connections"
    nodes = []
    edges = []
    index = 0
    for tup in global_vertex_mappings:
        nodes.append(vertices[tup[0]])
        nodes.append(vertices[tup[1]])
        edges.append([index, index + 1])
        index += 2

    nodes = np.asarray(nodes)
    edges = np.asarray(edges)

    ps.register_curve_network("virtual connections", nodes, edges, radius=0.001)

    ps.show()

if __name__ == "__main__":
    main()
