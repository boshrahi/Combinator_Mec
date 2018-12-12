import cc.redberry.combinatorics.Combinatorics;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class Main {

    public static void main(String[] args) throws IOException {
        final int numOfUsers = 100;
        final int numOfApps = 2; // just two for now
        final String graph = Graph.NOEL;
        //TODO make 2 python program to calculate all of answers and use it as a file

        // test
        for (int numOfVRCPerApp = 1 ; numOfVRCPerApp <= 4 ; numOfVRCPerApp++){
            Simulation2 simulation = new Simulation2(Graph.TEST, numOfVRCPerApp, numOfApps);  // just works for 2 apps for now
            List<String> all  = simulation.combination();
            Path file = Paths.get("test_vm"+numOfVRCPerApp+".txt");
            Files.write(file, all, Charset.forName("UTF-8"));
        }
        //  noel
        for (int numOfVRCPerApp = 1 ; numOfVRCPerApp <= 4 ; numOfVRCPerApp++){
            Simulation2 simulation = new Simulation2(Graph.NOEL, numOfVRCPerApp, numOfApps);  // just works for 2 apps for now
            List<String> all  = simulation.combination();
            Path file = Paths.get("noel_vm"+numOfVRCPerApp+".txt");
            Files.write(file, all, Charset.forName("UTF-8"));
        }
        //  shentel
        for (int numOfVRCPerApp = 1 ; numOfVRCPerApp <= 4 ; numOfVRCPerApp++){
            Simulation2 simulation = new Simulation2(Graph.SHENTEL, numOfVRCPerApp, numOfApps);  // just works for 2 apps for now
            List<String> all  = simulation.combination();
            Path file = Paths.get("shentel_vm"+numOfVRCPerApp+".txt");
            Files.write(file, all, Charset.forName("UTF-8"));
        }
        //  sago
        for (int numOfVRCPerApp = 1 ; numOfVRCPerApp <= 4 ; numOfVRCPerApp++){
            Simulation2 simulation = new Simulation2(Graph.SAGO, numOfVRCPerApp, numOfApps);  // just works for 2 apps for now
            List<String> all  = simulation.combination();
            Path file = Paths.get("sago_vm"+numOfVRCPerApp+".txt");
            Files.write(file, all, Charset.forName("UTF-8"));
        }
        //  spiralight
        for (int numOfVRCPerApp = 1 ; numOfVRCPerApp <= 4 ; numOfVRCPerApp++){
            Simulation2 simulation = new Simulation2(Graph.SPIRALIGHT, numOfVRCPerApp, numOfApps);  // just works for 2 apps for now
            List<String> all  = simulation.combination();
            Path file = Paths.get("spiralight_vm"+numOfVRCPerApp+".txt");
            Files.write(file, all, Charset.forName("UTF-8"));
        }
        //  missouri
        for (int numOfVRCPerApp = 1 ; numOfVRCPerApp <= 4 ; numOfVRCPerApp++){
            Simulation2 simulation = new Simulation2(Graph.MISSOURI, numOfVRCPerApp, numOfApps);  // just works for 2 apps for now
            List<String> all  = simulation.combination();
            Path file = Paths.get("missouri_vm"+numOfVRCPerApp+".txt");
            Files.write(file, all, Charset.forName("UTF-8"));
        }
    }


}
class Simulation2 {

    private String graphType;
    private int numOfVRCPerApp;
    private int numOfUsers;
    private int numOfApps;
    private GraphModel graph;
    private static List<VRCnodeModel> placementsForSingleApp;
    private static List<String> allVMS;
    //------------------------request


    Simulation2(String graphType, int numOfVRCPerApp, int numOfApps) {

        this.graphType = graphType;
        Graph graph = new Graph(graphType);
        this.graph = graph.getGraphModel();
        this.numOfApps = numOfApps;
        this.numOfVRCPerApp = numOfVRCPerApp;
    }

    /*
     * Algorithm 1 in paper
     * */
    List<String> combination() {

        List<String> all = initialAllPlacements(numOfVRCPerApp, numOfApps, graph);
        return all;
    }


    private List<String> initialAllPlacements(int numOfVRCPerApp, int numOfApps, GraphModel graph) {

        // System.out.println(enumeratePlacementCounts());
        //---------------------------------
        java.util.List<java.util.List> all = new ArrayList<>();

        for (int j = 0; j < numOfApps; j++) {

            placementsForSingleApp = new ArrayList<>();

            Combinatorics.combinations(graph.nodeNum, numOfVRCPerApp)
                    .stream()
                    .map(Simulation2::combinChangeIntToString)
                    .forEach(Simulation2::combinWritePlacementToList);
            all.add(placementsForSingleApp);


        }

        List<VRCnodeModel> list = all.get(0);
        VRCnodeModel model;
        VRCnodeModel[] set1 = new VRCnodeModel[list.size()];

        for (int i = 0; i < list.size(); i++) { // firstApp
            model = list.get(i);
            set1[i] = model;
        }
        List<VRCnodeModel> list2 = all.get(1);
        VRCnodeModel model2;
        VRCnodeModel[] set2 = new VRCnodeModel[list2.size()];
        for (int i = 0; i < list2.size(); i++) { // second App
            model2 = list2.get(i);
            set2[i] = model2;
        }
        allVMS = new ArrayList<>();
        Combinatorics.tuples(set1, set2)
                .stream()
                .map(Simulation2::tupleCombine)
                .forEach(Simulation2::tupleWriteToList);
        return allVMS;
    }

    private static void tupleWriteToList(String numbers) {
        allVMS.add(numbers);
    }

    private static String tupleCombine(VRCnodeModel[] vrCnodeModels) {
        String allVms = "";
        for (int i = 0; i < vrCnodeModels.length; i++) {
            VRCnodeModel model = vrCnodeModels[i];
            for (int j = 0; j < model.map.size(); j++) {
                long vm_placement = model.map.get(j);
                allVms = allVms + vm_placement + ",";

            }

        }
        String[] arr = allVms.split(",");
        return allVms;
    }

    private static String combinChangeIntToString(int[] ints) {
        String newString = "";
        for (int index = 0; index < ints.length; index++) {
            if (index == ints.length - 1) newString = newString + ints[index];
            else newString = newString + ints[index] + " ";
        }
        return newString;
    }

    private static void combinWritePlacementToList(String s) {
        String[] array = s.split(" ");
        VRCnodeModel placement = new VRCnodeModel();
        for (int vrcIndex = 0; vrcIndex < array.length; vrcIndex++) {
            placement.map.put(vrcIndex, Long.valueOf(array[vrcIndex]));
        }
        placementsForSingleApp.add(placement);

    }

    private long enumeratePlacementCounts() {
        long nodes = graph.nodeNum;
        long multi = 1;
        for (int i = 1; i <= numOfVRCPerApp; i++) {
            multi = multi * nodes;
            nodes--;
        } // enumerate places form VRCs in utils.Graph
        multi = (long) Math.pow(multi, numOfApps);
        return multi;
    }

}


/**
 * @author boshra
 * this class is for making graph of mec server: we use graphs of real world obtained from http://www.topology-zoo.org/explore.html
 * utils.Graph 1 : Spiralight               model.NodeModel = 15 Link = 16
 * utils.Graph 2 : Sago                     model.NodeModel = 18 Link = 17
 * utils.Graph 3 : Shentel                  model.NodeModel = 28 Link = 35
 * utils.Graph 4 : Missouri                 model.NodeModel = 67 Link = 83
 */

 class Graph {
    public static final String SPIRALIGHT = "SPIRALIGHT";
    public static final String SAGO = "SAGO";
    public static final String SHENTEL = "SHENTEL";
    public static final String MISSOURI = "MISSOURI";
    public static final String TEST = "TEST";
    public static final String NOEL = "NOEL";
    public GraphModel model;


    public Graph(String type) {
        model = new GraphModel();
        switch (type) {
            case SPIRALIGHT:
                model = creatGraph(SPIRALIGHT);
                break;
            case NOEL:
                model = creatGraph(NOEL);
                break;
            case SAGO:
                model = creatGraph(SAGO);
                break;
            case SHENTEL:
                model = creatGraph(SHENTEL);
                break;
            case MISSOURI:
                model = creatGraph(MISSOURI);
                break;
            case TEST:
                model = creatGraph(TEST);
                break;
        }


    }

    public GraphModel getGraphModel() {
        return model;
    }

    private GraphModel creatGraph(String type) {
        JSONParser parser = new JSONParser();
        JSONObject jsonObject = null;
        try {

            switch (type) {
                case TEST:
                    jsonObject = (JSONObject) parser.parse(new FileReader("test.json"));
                    break;
                case NOEL:
                    jsonObject = (JSONObject) parser.parse(new FileReader("Noel.json"));
                    break;
                case SPIRALIGHT:
                    jsonObject = (JSONObject) parser.parse(new FileReader("Spiralight.json"));
                    break;
                case SAGO:
                    jsonObject = (JSONObject) parser.parse(new FileReader("Sago.json"));
                    break;
                case SHENTEL:
                    jsonObject = (JSONObject) parser.parse(new FileReader("Shentel.json"));
                    break;
                case MISSOURI:
                    jsonObject = (JSONObject) parser.parse(new FileReader("Missouri.json"));
                    break;
            }
        } catch (Exception e) {
        }
        GraphModel graphModel = parseJsonFile(jsonObject);
        return graphModel;
    }

    private GraphModel parseJsonFile(JSONObject jsonObject) {

        GraphModel graphModel = new GraphModel();
        JSONArray nodesArr = (JSONArray) jsonObject.get("nodes");
        graphModel.nodeNum = nodesArr.size();

        JSONArray edgeArr = (JSONArray) jsonObject.get("edges");
        graphModel.linkNum = edgeArr.size();

        for (Object aNodesArr : nodesArr) {
            JSONObject object = (JSONObject) aNodesArr;
            String name = (String) object.get("label");
            long id = (long) object.get("id");
            NodeModel nodeModel = new NodeModel(name, id);
            graphModel.nodeModelList.add(nodeModel);

        }

        graphModel.prepareAdjacencyList();
        for (int i = 0; i < graphModel.nodeNum; i++) {
            for (int j = 0; j < graphModel.linkNum; j++) {
                JSONObject obj = (JSONObject) edgeArr.get(j);
                long source = (long) obj.get("source");
                long target = (long) obj.get("target");
                String id = (String) obj.get("id");
                long distance = (long) obj.get("distance");

                if (source == i) {
                    EdgeModel edgeModel = new EdgeModel(source, target, id, distance);
                    graphModel.addAdjacencyEdge((int) source, (int) target);
                    graphModel.edgeModelList.add(edgeModel);
                }

            }

        }

        return graphModel;
    }
}

 class GraphModel {
    public int nodeNum;
    public int linkNum;
    public List<EdgeModel> edgeModelList = new ArrayList<>();
    public List<NodeModel> nodeModelList = new ArrayList<>();
    public int maxLevel;

    private LinkedList<Integer> adjacency[]; //Adjacency Lists

    // Function to add an edge into the graph
    public void prepareAdjacencyList() {
        adjacency = new LinkedList[nodeNum];
        for (int i = 0; i < nodeNum; ++i) {
            adjacency[i] = new LinkedList();
        }
    }

    public void addAdjacencyEdge(int source, int target) {
        adjacency[source].add(target);
        adjacency[target].add(source);
    }

    public long[][] makeGraphMatrix() {
        long graph[][] = new long[nodeNum][nodeNum];
        for (int xIndex = 0; xIndex < nodeNum; xIndex++) {
            for (int yIndex = 0; yIndex < nodeNum; yIndex++) {
                for (int edge = 0; edge < edgeModelList.size(); edge++) {
                    if (xIndex == (int) edgeModelList.get(edge).source && yIndex == (int) edgeModelList.get(edge).target ||
                            xIndex == (int) edgeModelList.get(edge).target && yIndex == (int) edgeModelList.get(edge).source) {
                        graph[xIndex][yIndex] = edgeModelList.get(edge).distance;
                        break;
                    } else {
                        graph[xIndex][yIndex] = 0;
                    }
                }

            }
        }
        return graph;
    }

    // prints getNhops traversal from a given source s
    public HashMap<Integer, String> getNhops(int s) {
        // Mark all the vertices as not visited(By default
        // set as false)
        boolean visited[] = new boolean[nodeNum];

        // Create a queue for getNhops
        LinkedList<Integer> queue = new LinkedList<Integer>();

        HashMap<Integer, String> levels = new HashMap<>();
        for (int level = 0; level < nodeNum; level++) {
            levels.put(level, "");
        }

        // Mark the current node as visited and enqueue it
        visited[s] = true;
        queue.add(s);
        int level = 0;
        int maxLevel = 1;
        ArrayList<Integer> childList = new ArrayList<Integer>();
        while (queue.size() != 0) {
            // Dequeue a vertex from queue and print it
            s = queue.poll();
            if (!childList.contains(level)) level++;
            if (!childList.isEmpty()) childList.remove(0);

            // Get all adjacent vertices of the dequeued vertex s
            // If a adjacent has not been visited, then mark it
            // visited and enqueue it
            Iterator<Integer> i = adjacency[s].listIterator();
            while (i.hasNext()) {
                int n = i.next();
                if (!visited[n]) {
                    visited[n] = true;
                    queue.add(n);
                    childList.add(level + 1);
                    String hop_neighbours = levels.get(level);
                    if (!hop_neighbours.equals(""))
                        levels.put(level, hop_neighbours + "," + n);
                    else levels.put(level, String.valueOf(n));
                    maxLevel = level;

                }
            }
        }
        this.maxLevel = maxLevel;
        return levels;
    }


    // Function that implements Dijkstra's
    // single source shortest path
    // algorithm for a graph represented
    // using adjacency matrix
    // representation
    public ShortestPath dijkstra(long[][] adjacencyMatrix,
                                 int startVertex) {
        int nVertices = adjacencyMatrix[0].length;

        // shortestDistances[i] will hold the
        // shortest distance from src to i
        long[] shortestDistances = new long[nVertices];

        // added[i] will true if vertex i is
        // included / in shortest path tree
        // or shortest distance from src to
        // i is finalized
        boolean[] added = new boolean[nVertices];

        // Initialize all distances as
        // INFINITE and added[] as false
        for (int vertexIndex = 0; vertexIndex < nVertices;
             vertexIndex++) {
            shortestDistances[vertexIndex] = Integer.MAX_VALUE;
            added[vertexIndex] = false;
        }

        // Distance of source vertex from
        // itself is always 0
        shortestDistances[startVertex] = 0;

        // Parent array to store shortest
        // path tree
        long[] parents = new long[nVertices];

        // The starting vertex does not
        // have a parent
        parents[startVertex] = -1;

        // Find shortest path for all
        // vertices
        for (int i = 1; i < nVertices; i++) {

            // Pick the minimum distance vertex
            // from the set of vertices not yet
            // processed. nearestVertex is
            // always equal to startNode in
            // first iteration.
            int nearestVertex = -1;
            long shortestDistance = Integer.MAX_VALUE;
            for (int vertexIndex = 0;
                 vertexIndex < nVertices;
                 vertexIndex++) {
                if (!added[vertexIndex] &&
                        shortestDistances[vertexIndex] <
                                shortestDistance) {
                    nearestVertex = vertexIndex;
                    shortestDistance = shortestDistances[vertexIndex];
                }
            }

            // Mark the picked vertex as
            // processed
            added[nearestVertex] = true;

            // Update dist value of the
            // adjacent vertices of the
            // picked vertex.
            for (int vertexIndex = 0;
                 vertexIndex < nVertices;
                 vertexIndex++) {
                long edgeDistance = adjacencyMatrix[nearestVertex][vertexIndex];

                if (edgeDistance > 0
                        && ((shortestDistance + edgeDistance) <
                        shortestDistances[vertexIndex])) {
                    parents[vertexIndex] = nearestVertex;
                    shortestDistances[vertexIndex] = shortestDistance +
                            edgeDistance;
                }
            }
        }

        ShortestPath shortestPath = new ShortestPath();
        shortestPath.shorestDist = shortestDistances;
        shortestPath.parent = parents;
        return shortestPath;
        // printSolution(startVertex, shortestDistances, parents);
    }

    // A utility function to print
    // the constructed distances
    // array and shortest paths
    private void printSolution(int startVertex,
                               long[] distances,
                               long[] parents) {
        int nVertices = distances.length;
        System.out.print("Vertex\t Distance\tPath");

        for (int vertexIndex = 0;
             vertexIndex < nVertices;
             vertexIndex++) {
            if (vertexIndex != startVertex) {
                System.out.print("\n" + startVertex + " -> ");
                System.out.print(vertexIndex + " \t\t ");
                System.out.print(distances[vertexIndex] + "\t\t");
                printPath(vertexIndex, parents);
            }
        }
    }

    // Function to print shortest path
    // from source to currentVertex
    // using parents array
    private void printPath(long currentVertex,
                           long[] parents) {

        // Base case : Source node has
        // been processed
        if (currentVertex == -1) {
            return;
        }
        printPath(parents[(int) currentVertex], parents);
        System.out.print(currentVertex + " ");
    }
}
 class ShortestPath {
    public long shorestDist[];
    public long parent[];
}

class EdgeModel {

    public long source;
    public long target;
    public String id;
    public long distance;
    public EdgeModel(long source, long target, String id,long distance) {
        this.source = source;
        this.target= target;
        this.id = id;
        this.distance = distance;
    }
}
class NodeModel {
    public String name;
    public long id;
    public NodeModel(String name, long id){
        this.name = name;
        this.id = id;
    }
}
class VRCnodeModel { // every VRC mapped to a single model.NodeModel

    public HashMap<Integer , Long > map = new HashMap<>();

}