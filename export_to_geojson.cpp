#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>


std::map<std::string, std::vector<std::vector<std::pair<int, int>>>> polygons;

void parse_polygons(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line;
    std::string current_key;
    std::vector<std::pair<int, int>> current_polygon;
    std::istringstream iss;
    bool in_polygon = false;
    //while (file >> line) {
    while (std::getline(file, line)) {
        if (line[0] == '\'') {
            if (!current_key.empty()) {
                polygons[current_key].push_back(current_polygon);
                current_polygon.clear();
            }
//            current_key = line.substr(1, line.find_last_of('\'') - 1);
            current_key = line.substr(1, line.find('\'', 2) - 1);
        } else if (line[0] == 'p') {
            if (in_polygon) {
                polygons[current_key].push_back(current_polygon);
                current_polygon.clear();
            }
            in_polygon = true;
//        } else if (in_polygon && line.find_first_not_of("0123456789 ") == std::string::npos) {
//            std::istringstream coords(line);
//            int x, y;
//            coords >> x >> y;
//           current_polygon.push_back({x, y});
//        }
        } else if (in_polygon) {
            iss.clear();
            iss.str(line);
            int x, y;
            if (iss >> x >> y) {
                current_polygon.push_back({x, y});
            }
        }
    }

    if (!current_key.empty() && !current_polygon.empty()) {
        polygons[current_key].push_back(current_polygon);
    }
}

void export_geojson(const std::string& path) {
    std::ofstream out(path);
    out << "{\n";
    out << "  \"type\": \"FeatureCollection\",\n";
    out << "  \"features\": [\n";

    bool first = true;
    for (const auto& pair : polygons) {
        if (!first) {
            out << ",\n";
        }
        out << "    {\n";
        out << "      \"type\": \"Feature\",\n";
        out << "      \"properties\": {\n";
        out << "        \"name\": \"" << pair.first << "\"\n";
        out << "      },\n";
        out << "      \"geometry\": {\n";
        out << "        \"type\": \"MultiPolygon\",\n";
        out << "        \"coordinates\": [\n";

        bool first_polygon = true;
        for (const auto& polygon : pair.second) {
            if (!first_polygon) {
                out << ",\n";
            }
            out << "          [\n";
            bool first_point = true;
            for (const auto& point : polygon) {
                if (!first_point) {
                    out << ",\n";
                }
                out << "            [" << point.first << ", " << point.second << "]";
                first_point = false;
            }
            out << "          ]\n";
            first_polygon = false;
        }

        out << "        ]\n";
        out << "      }\n";
        out << "    }";
        first = false;
    }

    out << "\n  ]\n";
    out << "}";
}

int main() {
    float start = clock();
    std::string file_path = "/workspace/projs/dth/proj_lib/L13450163/check/run_dicectory/0722/DTH_PERI_RevC.drc.results"; 
    parse_polygons(file_path);

    std::string output_path = "output.geojson";
    //export_geojson(output_path);

    std::cout << "GeoJSON has been exported to " << output_path << std::endl;

    // 输出polygons map
    for (const auto& pair : polygons) {
        std::cout << "Key: " << pair.first << std::endl;
        std::cout << "Key_Length " << pair.second.size() << std::endl;
    //    for (const auto& polygon : pair.second) {
    //        for (const auto& point : polygon) {
    //            std::cout << "  (" << point.first << ", " << point.second << ")" << std::endl;
    //        }
    //    }
    }
    std::cout << "Time " << clock() - start << std::endl;
    return 0;
}
