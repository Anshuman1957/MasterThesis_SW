import numpy as np
import math
from EventModel import Event

def P2Length(P1:np.array,P2:np.array,rounding:int=4):
    '''Computes the distance between the 2 input points using the built-in math.dist function'''
    Distance = math.dist(P1,P2)
    Distance = round(Distance,rounding)
    return Distance

def P3Area(P1:np.array,P2:np.array,P3:np.array,rounding:int=4):
    '''Computes the area of a triangle given the 3 vertex coordinates. Uses Heron's semi-perimeter formula'''
    L1 = math.dist(P1,P2)
    L2 = math.dist(P2,P3)
    L3 = math.dist(P1,P3)
    S = (L1+L2+L3)/2
    Area = round(np.sqrt(S*(S-L1)*(S-L2)*(S-L3)),rounding)

    return Area

def InfoFactorA(Alpha:list[float]=None,Time:float=None,EventList:list[Event]=None,EventStates:list=None,Estimate:bool=False,Measure:bool=None,rounding:int=4,DEBUG=False):
    '''Computes the first term of the information factor given alpha, time, event list, and whether the event values are to be estimated or measured.
    Estimate means to measure the probability of the event(s) occuring at a given time.
    Measure means to get the true value from the object at the time.'''
    if Measure is None:
        Measure = not Estimate 
    t = Time

    if isinstance(Alpha,(int,float)):
        Alpha = np.array([Alpha] * len(EventList),dtype=float)
    else:
        Alpha = np.array(Alpha,dtype=float)
    
    Denominator = 1
    Numerator = 0

    if Estimate == True:
        for index,event in enumerate(EventList):
            Num = Alpha[index] * 2 * np.absolute(0.5 - event.CheckEventTrigger(Var=t,Estimate=True))
            Numerator += Num
    elif Measure == True:
        for index,event in enumerate(EventList):
            EventState = float(EventStates[index])
            Num = Alpha[index] * 2 * np.absolute(0.5 - event.CheckEventTrigger(Var=t,Measure=True))
            Numerator += Num
    else: # Both false, consider event hold states
        for index,event in enumerate(EventList):
            Num = Alpha[index] * 2 * np.absolute(0.5 - event.CheckEventTrigger(Var=t))
            Numerator += Num
    InformationFactorA = round(Numerator/Denominator,rounding)
    return InformationFactorA


def InfoFactorB(Beta:float=0.5,P1:np.array=None,P2:np.array=None,P_MAX:float=None,MIN_FACTOR:float=None,rounding:int=4):
    '''Computes the second term of the information factor given beta, measured and expected positions as 3D coordinates, minimum information factor \
        and maximum distance'''
    deltaDistance = P2Length(P1,P2)
    InformationFactorB = Beta * np.e**(-1*(deltaDistance)*MIN_FACTOR/P_MAX)
    InformationFactorB = round(InformationFactorB,rounding)

    return InformationFactorB

