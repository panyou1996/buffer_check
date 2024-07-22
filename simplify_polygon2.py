import shapely.geometry as sg
import matplotlib.pyplot as plt

def simplify_polygon(polygon):
    """
    Simplifies a polygon with only vertical and horizontal line segments 
    into a rectangle.

    Args:
        polygon: A shapely.geometry.Polygon object.

    Returns:
        A list of coordinates representing the simplified rectangle.
    """

    # Create a list of vertices from the polygon's exterior coordinates
    vertices = list(polygon.exterior.coords)

    # Sort the vertices by y-coordinate
    sorted_vertices = sorted(vertices, key=lambda v: v[1])

    # Create lists for horizontal and vertical segments
    horizontal_segments = []
    vertical_segments = []

    # Find horizontal segments (same y-coordinate)
    for i in range(len(sorted_vertices) - 1):
        if sorted_vertices[i][1] == sorted_vertices[i + 1][1]:
            horizontal_segments.append((sorted_vertices[i], sorted_vertices[i + 1]))

    # Find vertical segments (same x-coordinate)
    for i in range(len(sorted_vertices)):
        for j in range(i + 1, len(sorted_vertices)):
            if sorted_vertices[i][0] == sorted_vertices[j][0]:
                vertical_segments.append((sorted_vertices[i], sorted_vertices[j]))

    # Find the rectangle with the largest area
    max_area = 0
    max_rectangle = None
    for horizontal_segment in horizontal_segments:
        for vertical_segment in vertical_segments:
            rectangle = sg.Polygon([horizontal_segment[0], horizontal_segment[1], vertical_segment[1], vertical_segment[0]])
            area = rectangle.area
            if area > max_area:
                max_area = area
                max_rectangle = rectangle

    # Return the coordinates of the simplified rectangle
    return list(max_rectangle.exterior.coords)

def plot_polygon(vertices, color='blue'):
    """
    Plots a polygon using its exterior coordinates.

    Args:
        vertices: A shapely.geometry.Polygon object or a list of coordinates.
        color: The color of the plot (default: blue).
    """

    # Access exterior coordinates if vertices is a Polygon
    if isinstance(vertices, sg.Polygon):
        x = [v[0] for v in vertices.exterior.coords]
        y = [v[1] for v in vertices.exterior.coords]
    else:
        # Handle list of coordinates directly
        x = [v[0] for v in vertices]
        y = [v[1] for v in vertices]

    # Close the polygon for plotting
    x.append(x[0])
    y.append(y[0])
    plt.plot(x, y, color=color)

def plot_segments(segments, color='red'):
    """
    Plots line segments with markers.

    Args:
        segments: A list of line segments represented by pairs of coordinates.
        color: The color of the plot (default: red).
    """

    # Handle both list of segments and single segment cases
    if isinstance(segments, list):
        for segment in segments:
            x = [segment[0][0], segment[1][0]]
            y = [segment[0][1], segment[1][1]]
            plt.plot(x, y, color=color, marker='o')
    else:
        x = [segments[0][0], segments[1][0]]
        y = [segments[0][1], segments[1][1]]
        plt.plot(x, y, color=color, marker='o')

# Example usage
polygon = sg.Polygon([(0,10),(5,10),(5,0),(5.2,0),(5.2,8),(10,8),(10,-2),(10.3,-2),(10.3,3),(10.3,24),(5,24),(5,23.5),(0,23.5),(0,10)])
simplified_polygon = simplify_polygon(polygon)
print(simplified_polygon)

# Plot the results
fig, axs = plt.subplots(1, 2, figsize=(10, 5))

# Input polygon
axs[0].set_title('Input Polygon')
plot_polygon(polygon, color='blue')
axs[0].plot(*polygon.exterior.xy, color='blue')

# Simplified polygon
axs[1].set_title('Simplified Polygon')
plot_polygon(simplified_polygon, color='red')
axs[1].plot(*zip(*simplified_polygon), color='red')

plt.show()
