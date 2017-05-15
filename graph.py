""" A Python Class
A simple Python graph class, demonstrating the essential 
facts and functionalities of graphs.
"""

class Graph(object):

    def __init__(self, graph_dict=None):
        """ initializes a graph object 
            If no dictionary or None is given, 
            an empty dictionary will be used
        """
        if graph_dict == None:
            graph_dict = {}
        self.__graph_dict = graph_dict

    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []

    def add_edge(self, edge):
        """ assumes that edge is of type set, tuple or list; 
            between two vertices can be multiple edges! 
        """
        # edge = set(edge)
        (vertex1, vertex2) = tuple(edge)
        if vertex1 in self.__graph_dict:
            self.__graph_dict[vertex1].append(vertex2)
        else:
            self.__graph_dict[vertex1] = [vertex2]

    def __generate_edges(self):
        """ A static method generating the edges of the 
            graph "graph". Edges are represented as sets 
            with one (a loop back to the vertex) or two 
            vertices 
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                if {neighbour, vertex} not in edges:
                    edges.append({vertex, neighbour})
        return edges

    def __str__(self):
        res = "vertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res

    def breadth_first_search(self, root):
        rotate_order = []
        # TODO Fazer isso generico para o numero
        parent = []
        visited_node = []
        for x in xrange(0,20):
            parent.append(-1) 
            visited_node.append(0)
        visited_node[root] = 1
        queue = []
        queue.append(root)
        while len(queue) != 0 :
            v = queue[0]
            for neighbour in self.__graph_dict[v]:
                # print(queue)
                if visited_node[neighbour] == 0:
                    visited_node[neighbour] = 1
                    queue.append(neighbour)
                    rotate_order.append(neighbour)
                    parent[neighbour] = v
                elif neighbour in queue:
                    pass
                    # print('estou aqui')
                    #visitar aresta
            del queue[queue.index(v)]        
        return { 'order':rotate_order, 'parent':parent}

    def find_path(self, start_vertex, end_vertex, path=None):
        """ find a path from start_vertex to end_vertex 
            in graph """
        if path == None:
            path = []
        graph = self.__graph_dict
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            return path
        if start_vertex not in graph:
            return None
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_path = self.find_path(vertex, 
                                               end_vertex, 
                                               path)
                if extended_path: 
                    return extended_path
        return None

    def check_adjacent_face_has_vertices(self, faces, index, v1,v2, polygons):
        i = 0;
        l = []
        for face in faces:
            if i != index:
                count = 0;
                l = []
                for vertice in face:
                    if vertice == v1 or vertice == v2:
                        count = count + 1
                        l.append(vertice)
                if count == 2:                    
                    polygons[index].edges[i] = l
                    return i;
            i = i+1;

    def generate_graph(self, g, edges, polygons):
        i = 0;
        faces = [];
        for edge in edges :
            faces.append(edges[i][0].tolist())
            g[i] = [] 
            i = i+1 
        
        i=0;
        for face in faces:
            for vertice in range(len(face)-1):
                viz = self.check_adjacent_face_has_vertices(faces,i,face[vertice],face[vertice+1],polygons)
                self.add_edge([i,viz])
            viz = self.check_adjacent_face_has_vertices(faces,i,face[vertice+1],face[0],polygons)
            self.add_edge([i,viz])   
            i=i+1
        #print(g)    
        return

if __name__ == "__main__":

    g = { "a" : ["d", "f"],
      "b" : ["c"],
      "c" : ["b", "c", "d", "e"],
      "d" : ["a", "c"],
      "e" : ["c"],
      "f" : ["d"]
    }

    graph = Graph(g)

    print("Vertices of graph:")
    print(graph.vertices())

    print("Edges of graph:")
    print(graph.edges())

    print("Add vertex:")
    graph.add_vertex("z")

    print("Vertices of graph:")
    print(graph.vertices())
 
    print("Add an edge:")
    graph.add_edge({"a","z"})
    
    print("Vertices of graph:")
    print(graph.vertices())

    print("Edges of graph:")
    print(graph.edges())

    print('Adding an edge {"x","y"} with new vertices:')
    graph.add_edge({"x","y"})
    print("Vertices of graph:")
    print(graph.vertices())
    print("Edges of graph:")
    print(graph.edges())
    print(g)