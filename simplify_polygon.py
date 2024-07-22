import random
import matplotlib.pyplot as plt

def generate_polygon_with_holes(grid_size=10, num_notches=5):
    vertices = []
    
    # Start at a random position
    x, y = random.randint(0, grid_size), random.randint(0, grid_size)
    vertices.append((x, y))
    
    # Create outer boundary
    for _ in range(grid_size * 2):
        direction = random.choice(['horizontal', 'vertical'])
        step = random.randint(1, 3)
        if direction == 'horizontal':
            x = (x + step) % grid_size
        else:
            y = (y + step) % grid_size
        vertices.append((x, y))

    # Create notches
    for _ in range(num_notches):
        notch_x, notch_y = random.randint(1, grid_size - 1), random.randint(1, grid_size - 1)
        vertices.append((notch_x, notch_y))
        vertices.append((notch_x, notch_y + random.randint(1, 3)))
        vertices.append((notch_x + random.randint(1, 3), notch_y + random.randint(1, 3)))
        vertices.append((notch_x + random.randint(1, 3), notch_y))

    # Ensure the polygon is closed
    vertices.append(vertices[0])
    
    return vertices

def simplify_polygon(vertices):
    simplified_vertices = [vertices[0]]
    
    for i in range(1, len(vertices) - 1):
        prev_vertex = vertices[i - 1]
        current_vertex = vertices[i]
        next_vertex = vertices[i + 1]
        
        if (prev_vertex[0] == current_vertex[0] == next_vertex[0]) or (prev_vertex[1] == current_vertex[1] == next_vertex[1]):
            continue
        else:
            simplified_vertices.append(current_vertex)
    
    simplified_vertices.append(vertices[-1])
    
    return simplified_vertices

def vertices_to_segments(vertices):
    segments = []
    for i in range(len(vertices) - 1):
        segments.append((vertices[i], vertices[i + 1]))
    return segments

def plot_polygon(vertices, color='blue'):
    x = [v[0] for v in vertices]
    y = [v[1] for v in vertices]
    x.append(vertices[0][0])
    y.append(vertices[0][1])
    plt.plot(x, y, color=color)

def plot_segments(segments, color='red'):
    for segment in segments:
        x = [segment[0][0], segment[1][0]]
        y = [segment[0][1], segment[1][1]]
        plt.plot(x, y, color=color, marker='o')

# 生成复杂的多边形
polygon_vertices = generate_polygon_with_holes(grid_size=10, num_notches=5)
polygon_vertices = [(0,10),(5,10),(5,0),(5.2,0),(5.2,8),(10,8),(10,-2),(10.3,-2),(10.3,3),(10.3,24),(5,24),(5,23.5),(0,23.5),(0,10)]
simplified_vertices = simplify_polygon(polygon_vertices)
segments = vertices_to_segments(simplified_vertices)

print("Simplified Vertices:")
for vertex in simplified_vertices:
    print(vertex)

print("\nSegments:")
for segment in segments:
    print(segment)

# 作图展示
plt.figure(figsize=(10, 5))

# 输入的polygon
plt.subplot(1, 2, 1)
plot_polygon(polygon_vertices, color='blue')
plt.title("Input Polygon")
plt.gca().set_aspect('equal', adjustable='box')

# 输出的线段
plt.subplot(1, 2, 2)
plot_segments(segments, color='red')
plt.title("Simplified Segments")
plt.gca().set_aspect('equal', adjustable='box')

plt.show()
