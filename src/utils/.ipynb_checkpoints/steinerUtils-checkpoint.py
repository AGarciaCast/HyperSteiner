
from src.utils.fullSteinerSolverHyperbolic import steinerPoint3Hyp, steinerPoints4Hyp, DISTANCE_HYP, hyperbolicInnerAngleTriangle
from src.utils.fullSteinerSolverEuclidean import steinerPoint3Euc, steinerPoints4Euc, l2Distance, euclideanInnerAngleTriangle
from math import inf #used for Dijkstra functions

import numpy as np

DISTANCE_F = {**DISTANCE_HYP, 
             "Euclidean": l2Distance}


def innerAngleTriangle(u, v, w, space="Klein"):
    if space not in ["Klein", "Euclidean"]:
        raise ValueError("space should be either 'Klein' or 'Euclidean'")
    
    if space=="Klein":
        result = hyperbolicInnerAngleTriangle(u, v, w, model="Klein")
    else:
        result = euclideanInnerAngleTriangle(u, v, w)
    
    return result


def steinerPoint3(vert, model="Klein", dist2Points=1e-1, precise=True):

    if model in ["Klein", "Half"]:
        result = steinerPoint3Hyp(vert, model, precise=precise, dist2Points=dist2Points)
    elif model == "Euclidean":
        result = steinerPoint3Euc(vert, dist2Points=dist2Points)
    else:
        raise ValueError("model should be either 'Klein', 'Half', or 'Euclidean'")
    
    return result


def steinerPoints4(vert, topo, model="Klein", nIters=100, convDiff=1e-2, dist2Points=1e-1, precise=True):

    if model in ["Klein", "Half"]:
        result = steinerPoints4Hyp(vert, topo, model, nIters=nIters, convDiff=convDiff,
                                   precise=precise, dist2Points=dist2Points)
    elif model == "Euclidean":
        result = steinerPoints4Euc(vert, topo, nIters=nIters, convDiff=convDiff, dist2Points=dist2Points)
    else:
        raise ValueError("model should be either 'Klein', 'Half', or 'Euclidean'")
    
    return result



def steinerRatio(vert, mst, model="Klein", idxTerminals=None, idxSteinerPoint=0,
                 nIters=100, convDiff=1e-2, dist2Points = 1e-1, precise=True):
    
    if idxSteinerPoint is None:
        idxTerminals = [i for i in range(len(vert))]
        
    numVert = len(vert)
    if numVert==3:
        return steinerRatio3(vert=vert,
                             mst=mst,
                             idxTerminals=idxTerminals,
                             model=model,
                             idxSteinerPoint=idxSteinerPoint,
                             precise=precise,
                             dist2Points=dist2Points
                             )
    elif numVert==4:
        return steinerRatio4(vert=vert,
                             mst=mst,
                             idxTerminals=idxTerminals,
                             model=model,
                             idxSteinerPoint=idxSteinerPoint,
                             nIters=nIters, 
                             convDiff=convDiff,
                             precise=precise,
                             dist2Points=dist2Points
                             )

def steinerRatio3(vert, mst, idxTerminals, model="Klein", idxSteinerPoint=0, dist2Points=1e-1, precise=True):
    
    steinerPoint = steinerPoint3(vert, model, precise=precise, dist2Points=dist2Points)
    
    if steinerPoint is None:
        ratio = 1.0
        topology = None
        
    else:

        smt = 0
        for i in range(3):
            smt += DISTANCE_F[model](steinerPoint, vert[i])

        ratio = smt / mst
        
        if ratio >=1:
            ratio = 1.0
            steinerPoint = None
            topology = None
        else:
            topology = [[f"T{idxTerminals[i]}", f"S{idxSteinerPoint}"] for i in range(3)]

    return ratio, steinerPoint, topology



def steinerRatio4(vert, mst, idxTerminals, model="Klein", idxSteinerPoint=0, nIters=100,
                  convDiff=1e-2, dist2Points=1e-1, precise=True):
    topologiesIdx = [[[0,1], [2, 3]], [[0,2], [1, 3]]]
    bestSteinerPoints = None
    bestSMT = inf
    bestTopo = None
    
    for topo in topologiesIdx:
        steinerPoints = steinerPoints4(vert, topo, model=model,
                                        nIters=nIters,
                                        convDiff=convDiff,
                                        dist2Points=dist2Points,
                                        precise=precise)
        
        
        if steinerPoints is None:
            ratio = 1.0
            fstTopology = None
        
        else:

            smt = DISTANCE_F[model](steinerPoints[0], steinerPoints[1])
            for i in range(2):
                for j in range(2):
                    smt += DISTANCE_F[model](steinerPoints[i], vert[topo[i][j]])
                        
            if smt < bestSMT:
                bestSMT = smt 
                bestSteinerPoints = steinerPoints
                bestTopo = topo
            
            
        if bestSteinerPoints is not None:  

            ratio = bestSMT / mst
            
            if ratio >=1:
                ratio = 1.0
                bestSteinerPoints = None
                fstTopology = None
            else:
                
                fstTopology = [[f"T{idxTerminals[bestTopo[i][j]]}", f"S{idxSteinerPoint + i}"] for i in range(2) for j in range(2)]
                fstTopology.append([f"S{idxSteinerPoint}", f"S{idxSteinerPoint + 1}"])
    
    
    return ratio, bestSteinerPoints, fstTopology

