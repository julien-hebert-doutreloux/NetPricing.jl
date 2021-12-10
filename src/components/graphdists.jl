dists_from(graph, o) = dijkstra_shortest_paths(graph, o).dists
dists_to(graph, d) = dijkstra_shortest_paths(reverse(graph), d).dists
