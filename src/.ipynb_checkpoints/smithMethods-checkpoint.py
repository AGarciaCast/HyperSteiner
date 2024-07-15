from src.utils.steinerUtils import *
from src.utils.fullSteinerSolverHyperbolic import *
from src.utils.delaunay2d import compute_Voronoi_Delaunay
from src.utils.graphsUtils import DisjSet
from queue import PriorityQueue
from collections import defaultdict 




def buildQueue(points, triangulationListEdges, triangulationListTriang, mapEdge2Triangle,
               space="Klein", maxgroup=4, nIters=100, convDiff=1e-2, dist2Points=1e-1, precise=True):
    
    verticesDict = dict()
    for i, p in enumerate(points):
        verticesDict[f"T{i}"] = p
    
    numPoints = len(points)
    # Building the edges-distance map:
    edgesWeight = []
    for (edge_idx0, edge_idx1) in triangulationListEdges:
            p0 = points[edge_idx0]
            p1 = points[edge_idx1]
            dist = DISTANCE_F[space](p1, p0)
            edgesWeight.append([f"T{edge_idx0}", f"T{edge_idx1}", dist])

    edgesWeight = sorted(edgesWeight, 
                        key=lambda item: item[2]) 
    
    
    auxGraphDisjSet = DisjSet(numPoints)
    # This will store the resultant MST 
    mstGraph = [] 
    adjTriangInfo = {"numEdges": np.zeros([len(triangulationListTriang)]),
                     "mstW": np.zeros([len(triangulationListTriang)]),
                     "dictEdges": defaultdict(list)
    }
    orderedSelectedSubsets = PriorityQueue()
    selectedTriangles = []

    # An index variable, used for sorted edges 
    i = 0

    # An index variable, used for result[] 
    e = 0
    
    idxSteinerPoint=0
    fullMSTval = 0
    # Number of edges to be taken is less than to V-1 
    while e < numPoints - 1: 

        # Pick the smallest edge and increment 
        # the index for next iteration 
        u, v, w = edgesWeight[i] 
        i = i + 1
        x = auxGraphDisjSet.find(u) 
        y = auxGraphDisjSet.find(v) 

        # If including this edge doesn't 
        # cause cycle, then include it in result 
        # and increment the index of result 
        # for next edge 
        if x != y: 
            e = e + 1
            fullMSTval += w
            mstGraph.append([[u, v]]) 
            u = int(u[1:])
            v = int(v[1:])
            try:
                adjTraingles = mapEdge2Triangle[(u, v)]
            except KeyError:
                adjTraingles = mapEdge2Triangle[(v, u)]
            
                
            for tr in adjTraingles:
                idxPointsTriangle = triangulationListTriang[tr]
                # Edges with no triangles are skiped
                if len(idxPointsTriangle)>2:
                    adjTriangInfo["numEdges"][tr] += 1
                    adjTriangInfo["mstW"][tr] += w
                    adjTriangInfo["dictEdges"][tr].append((u, v, w))
                    
                    if adjTriangInfo["numEdges"][tr] == 2:
                        
                        ratio, steinerPoint, topology  = steinerRatio(vert=points[idxPointsTriangle],
                                                            mst=adjTriangInfo["mstW"][tr],
                                                            model=space,
                                                            idxTerminals=idxPointsTriangle,
                                                            idxSteinerPoint=idxSteinerPoint,
                                                            dist2Points=dist2Points,
                                                            precise=precise)
                        
                      
                        
                        
                        if steinerPoint is not None and not np.isnan(ratio):
                            verticesDict[f"S{idxSteinerPoint}"] = steinerPoint
                            idxSteinerPoint += 1
                            orderedSelectedSubsets.put((ratio, topology))
                            selectedTriangles.append(tr)
            
            auxGraphDisjSet.union(x, y) 
        # Else discard the edge 
    
    if maxgroup==4:
        for tr in selectedTriangles:
            idxPointsTriangle = triangulationListTriang[tr]
            for i in range(3):
                            
                u = idxPointsTriangle[i]
                v = idxPointsTriangle[(i+1)%3]
                w = idxPointsTriangle[(i+2)%3]
                
                try:
                    trAdj = list(set(mapEdge2Triangle[(u, v)]) - {tr})
                except KeyError:
                    trAdj = list(set(mapEdge2Triangle[(v, u)]) - {tr})
               
                # Not boundary edge
                if len(trAdj) > 0:
                    trAdj = trAdj[0]
                    
                    # check if is a frond
                    if adjTriangInfo["numEdges"][tr] + adjTriangInfo["numEdges"][trAdj] > 2: 
                        verticesFrond = [w, u, v] + list(set(triangulationListTriang[trAdj])-{u,v}) 
                        mstFrond = sum(e[2] for e in set(adjTriangInfo["dictEdges"][tr] + adjTriangInfo["dictEdges"][trAdj]))
                        ratio, steinerPoints, topology  = steinerRatio(vert=points[verticesFrond],
                                                        mst=mstFrond,
                                                        idxTerminals=verticesFrond,
                                                        idxSteinerPoint=idxSteinerPoint,
                                                        model=space,
                                                        nIters=nIters,
                                                        convDiff=convDiff,
                                                        dist2Points=dist2Points,
                                                        precise=precise)
                        
                        if steinerPoints is not None and not np.isnan(ratio):
                            for p in steinerPoints:
                                verticesDict[f"S{idxSteinerPoint}"] = p
                                idxSteinerPoint += 1
                                
                            orderedSelectedSubsets.put((ratio, topology))
                        
        
    auxCons = []
    
    while not orderedSelectedSubsets.empty():
        auxCons.append(orderedSelectedSubsets.get()[1])
    

    constructionList = auxCons + mstGraph
    mstGraph = [el[0] for el in mstGraph]
    
    return constructionList, verticesDict, fullMSTval, mstGraph


def greedyConcatenation(points, constructionList, verticesDict, space="Klein"):
    # An index variable, used for sorted set of points
    i = 0

    vertex2Dict = dict()
    visitedVertex = defaultdict(lambda: False)
    smt = 0
    auxGraphDisjSet = DisjSet(len(points))
    resultGraph = []
    renameDict=dict()
    counterSteiner = 0
    numFST3 = 0
    numFST4 = 0
    
    while auxGraphDisjSet.numConnectedComponents() > 1 and i < len(constructionList):  

        newEdges = constructionList[i] 
        i = i + 1
        
        fathers = []
        auxSet = set()
     
        for e in newEdges:
            u = e[0]
            v = e[1]
            if u not in auxSet:
                fathers.append(auxGraphDisjSet.find(u))
                auxSet.add(u)
                
            if v not in auxSet:
                fathers.append(auxGraphDisjSet.find(v)) 
                auxSet.add(v)

        
        if  len(fathers) == len(set(fathers)):
            numEdges = len(newEdges)
            numFST3 += int(numEdges==3)
            numFST4 += int(numEdges==5)
            
            for e in newEdges:
                
                u = e[0]
                x = auxGraphDisjSet.find(u)
                
                if not visitedVertex[u]:
                    if u[0]=="S":
                        renameDict[u] = f"S{counterSteiner}"
                        counterSteiner += 1
                    else:
                        renameDict[u] = u
                    
                    vertex2Dict[renameDict[u]] = verticesDict[u]
                    visitedVertex[u] = True
                
                v = e[1]
                y = auxGraphDisjSet.find(v) 
                
                if not visitedVertex[v]:
                    if v[0]=="S":
                        renameDict[v] = f"S{counterSteiner}"
                        counterSteiner += 1
                    else:
                        renameDict[v] = v
                        
                    vertex2Dict[renameDict[v]] = verticesDict[v]
                    visitedVertex[v] = True
                
                smt += DISTANCE_F[space](verticesDict[u], verticesDict[v])
                
                
                resultGraph.append([renameDict[u], renameDict[v]])
                auxGraphDisjSet.union(x, y) 
   
    
    return resultGraph, vertex2Dict, smt, numFST3, numFST4


def heuristicMSTpipeline(points, space="Klein", triangMeth="DT", maxgroup=4,
                         nIters=100, convDiff=1e-2, dist2Points=1e-1, precise=True):
    """_summary_

    Args:
        points (_type_): _description_
        space (str, optional): "Klein" or "Euclidean". Defaults to "Klein".
        triangMeth (str, optional): _description_. Defaults to "DT".
        maxgroup (int, optional): _description_. Defaults to 4.
    """
    N = len(points)
    if N<3:
        resultGraph=list()
        if N==2:
            resultGraph.append(["T0", "T1"])
        
        verticesDict = dict()
        for i in range(N):
            verticesDict[f"T{i}"] = points[i]
            
        
        return resultGraph, verticesDict, 1, None, None
    
    if space not in ["Klein", "Euclidean"]:
        raise ValueError("space should be either 'Klein' or 'Euclidean'")
    
    if triangMeth not in ["DT", "GG"]:
        raise ValueError("space should be either 'DT' or 'GG'")
    
    if maxgroup not in [3, 4]:
        raise ValueError("maxgroup should be 3 or 4")

    if triangMeth == "DT":
        #{'original_points': S, 'Voronoi_points': V, 'Delaunay_triangulation': tri_list, 'Delaunay_edges': edges,
        # 'map_DelaunayEdge2Triangle': edge_map, 'map_DataPoint2VoronoiCell': voronoi_cell_map} 
        computedVoronoiDelaunay = compute_Voronoi_Delaunay(points, space = space)
        triangulationData = {"mapEdge2Triangle": computedVoronoiDelaunay["map_DelaunayEdge2Triangle"],
                             "triangulationListTriang": computedVoronoiDelaunay["Delaunay_triangulation"],
                             "triangulationListEdges": computedVoronoiDelaunay["Delaunay_edges"]
        }
        
    else:
        #TODO: obtain gabriel graph
        raise ValueError("triangMeth should be DT")

    
    constructionList, verticesDict, fullMSTval, mstGraph = buildQueue(points=points,
                                                                    space=space,
                                                                    maxgroup=maxgroup,
                                                                    nIters=nIters,
                                                                    convDiff=convDiff,
                                                                    dist2Points=dist2Points,
                                                                    precise=precise,
                                                                    **triangulationData)

    
    resultGraph, verticesDict, steinerVal, numFST3, numFST4 = greedyConcatenation(points=points,
                                                                 constructionList=constructionList,
                                                                 verticesDict=verticesDict,
                                                                 space=space
                                                                 )
   
    
    
    return  resultGraph, verticesDict, steinerVal/fullMSTval, mstGraph, triangulationData["triangulationListEdges"], numFST3, numFST4
