import re

def parse_document(file_path):
    with open(file_path, 'r') as file:
        data = {}
        current_key = None
        current_points = []
        polygon_pattern = re.compile(r"'([^']*)'")
        point_pattern = re.compile(r'(-?\d+\.?\d+)\s+(-?\d+\.?\d*)$')

        for line in file:
            line = line.strip()
            polygon_match = polygon_pattern.search(line)
            if polygon_match:
                current_key = polygon_match.group(1)
                if current_key and current_points:
                    if current_key not in data:
                        data[current_key] = []
                    data[current_key].append(current_points)
                    current_points = []
            elif line == current_key:
                if current_points:
                    data[current_key].append(current_points)
                    current_points = []
            elif line.startswith('p '):
                if current_points:
                    if current_key not in data:
                        data[current_key] = []
                    data[current_key].append(current_points)
                    current_points = []
            else:
                point_match = point_pattern.match(line)
                if point_match:
                    x, y = map(float, point_match.groups())
                    current_points.append((x/10000, y/10000))

        if current_key and current_points:
            if current_key not in data:
                data[current_key] = []
            data[current_key].append(current_points)

        return data

def format_output(data):
    for key, polygons in data.items():
        print("{}={}".format(key, polygons))

def output_polygons(data):
    return data

if __name__ == "__main__":
    file_path = r'/workspace/projs/dth/proj_lib/L13450163/check/run_dicectory/buffer_skill.drc.results'
    data = parse_document(file_path)
    formatted_data = output_polygons(data)
    print("read polygons datas finished")
    format_output(formatted_data)
