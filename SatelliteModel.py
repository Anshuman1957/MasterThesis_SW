import numpy as np
import matplotlib.pyplot as plt
from OrbitalElements import Orbit,Rotation
import sys

class Satellite:

    def __init__(self,ID:int = 0,InitialPosition:np.array = np.array([0,0,0]),InitialRelativePosition:np.array = np.array([0,0,0])):
        '''Initializes the Satellite with the following information:
        1. Identification Number for the satellite (Default 0)
        2. Co-ordinate position (Defaults to 3D co-ordinate geometry. Can choose quaternion representation)
        3. Initializes the parameter values (TBD)'''
        self.ID = ID # Identifier number for the satellite in the swarm
        self.FoV = 0 #Field of View: Visual/Communication range of the satellite in 3 dimensional space in metres^3
        self.PointingAccuracy = 1 # arcsecond, relevant for optical communication
        self.DeltaV_MAX = 1.53 # m/s, reference: Operational_Concept_of_PicoSat_Release_From_LEOSat.pdf
        self.Mass = 1 # kg, Pico-satellite
        self.UplinkRateMax = 1000 # kbps, tentative
        self.DownlinkRateMax = 1000 # kbps, tentative
        self.TransmitPower = 1 # Watts, typically 1 to 5 W, reference: R-REP-SA.2312-2014-PDF-E.pdf
        self.InitializePosition(InitialPosition) # Initialize the position of the satellite
        self.InitializeRelativePosition(InitialRelativePosition) # Initialize the relative position of the satellite
        self.Velocity = np.array([0,0,0]) # Velocity of the satellite in (x,y,z) in m/s
        self.Events = [] # Initialize an empty array of events attached to this satellite
        self.EventStates = [] # Current state of events (boolean)
        self.PreviousEventCheckTime = None # Simulation time in minutes when the last event check occured
        self.SetTriggerCondition(0.5)
        #self.BatteryCapacity =
        #self.BatteryConsumptionUplink = 
        #self.BatteryConsumtionDownlink = 
        #self.CommunicationProtocol =
        #self.RotationMatrix =
        #self.TranslationMatrix = 
        #self.Type = # Satellite role in the heterogeneous swarm

    @property
    def Position(self):
        '''Property Position: Value: Absolute 3-D coordinates of the satellite'''
        return self._Coordinates
    
    @property
    def RelativePosition(self):
        '''Property Relative Position: Relative (to Leader) 3-D coordinates of the satellite'''
        return self._RelativeCoordinates
    
    @property
    def PositionHistory(self):
        '''Property Position History: Value: List of absolute 3-D coordinates where the satellite existed'''
        return self._PositionHistory
    
    @property
    def RelativePositionHistory(self):
        '''Property Relative Position History: Value: List of relative 3-D coordinates where the satellite existed'''
        return self._RelativePositionHistory
    
    def ClearPositionHistory(self):
        '''Clears the position history of the satellite'''
        self._PositionHistory = []
        
    def ClearRelativePositionHistory(self):
        '''Clears the relative position history of the satellite'''
        self._RelativePositionHistory = []
    
    @Position.setter
    def Position(self,NewPosition:np.array):
        '''Property Position: Setter: Sets position to input position and updates position history'''
        self._Coordinates = NewPosition
        self._PositionHistory.append(self._Coordinates)
        
    @RelativePosition.setter
    def RelativePosition(self,NewPositionRelative:np.array):
        '''Property Relative Position: Setter: Sets relative position to input position and updates relative position history'''
        self._RelativeCoordinates = NewPositionRelative
        self._RelativePositionHistory.append(self._RelativeCoordinates)
    
    @property
    def TriggerThreshold(self):
        '''Property TriggerThreshold: float: Event trigger threshold to trigger transmission'''
        return self._TriggerThreshold
    
    @TriggerThreshold.setter
    def TriggerThreshold(self,NewThreshold=0):
        '''Property TriggerThreshold: Setter: Sets trigger threshold to input NewThreshold value'''
        self._TriggerThreshold = NewThreshold
    
    def InitializePosition(self,InitPosition:np.array=np.array([0,0,0])):
        '''Initializes initial position of the satellite'''
        self.ClearPositionHistory()
        self.Position = InitPosition
        
    def InitializeRelativePosition(self,InitPosition:np.array=np.array([0,0,0])):
        '''Initializes initial position of the satellite'''
        self.ClearRelativePositionHistory()
        self.RelativePosition = InitPosition

    def PlotPath(self,duration=1000,dt=100,velocityFunction=lambda x:np.array([x,x,x]),SaveFigure=False):
        '''Computes and plots the path given the input velocity functions/vectors and time parameters.
         Supported keyword arguments:
         dt : Time step - Defaults to 10s
         duration : Total time - Defaults to 1000s
         velocityFunction : List of constant Velocity functiosn acting on the satellite.
         SaveFigure : Boolean to save the intermediate scatter plots. Default False.
         '''
        
        for timeStep in range(0,duration,dt):
            self.Position = self.Position + dt*(velocityFunction(timeStep))
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('X-Axis')
        ax.set_ylabel('Y-Axis')
        ax.set_zlabel('Z-Axis')
        ax.set_xlim([0,500000])
        ax.set_ylim([0,45000])
        ax.set_zlim([0,1000])
        counter = 1
        for point in self.PositionHistory:
            ax.scatter(point[0],point[1],point[2],marker='o',sizes=[16])
            if SaveFigure == True:
                plt.savefig(str(counter)+".jpg")
            counter+=1
        plt.show()
        return (ax,fig)
    
    def DefineOrbit(self,Orbit:Orbit):
        '''Defines the Orbit of the satellite'''
        self.Orbit = Orbit

    def AddEvents(self,events):
        '''Adds the input Event(s) to the satellite event handler. Sets the initial event state to False.'''
        if isinstance(events,list):
            for event in events:
                self.Events.append(event)
                self.EventStates.append(False)
        elif isinstance(events,Event):
            self.Events.append(events)
            self.EventStates.append(False)
        else:
            raise TypeError
        
        return None

    def ClearEvents(self):
        '''Clears all attached events to the satellite'''
        self.Events = []
        self.EventStates = []

        return None

    def CheckEventStates(self,SimulationTime):
        '''Check if any events trigger based on simulation time. Returns a list of boolean values corresponding to the event trigger state.'''
        #deltaT = SimulationTime - self.PreviousEventCheckTime if self.PreviousEventCheckTime is not None else (SimulationTime - 0)
        T = SimulationTime
        self.PreviousEventCheckTime = T
        retVal = []
        for event in self.Events:
            EventState = event.CheckEventTrigger(T)
            retVal.append(EventState)
        self.CurrentEventStates = retVal
        return retVal
    
    def ArbitrateEventStates(self,Events=None,EventStates=None):
        '''Compute any external action required based on the event states. Call TriggerTransmission() if event(s) exceed defined action threshold.'''
        if Events is None:
            Events = self.Events
        if EventStates is None:
            EventStates = self.CurrentEventStates

        Priorities = np.array([event.Priority for event in Events])
        EventStatesInt = np.array(EventStates,dtype=int)
        TriggerFactor = np.sum(Priorities*EventStatesInt)
        if self.CheckTransmissionConditions(TriggerFactor):
            self.SetTriggerCondition(999)

    
    def TriggerTransmission(self):
        '''Trigger a (short-range radio) transmission broadcast intended for the leader of the swarm based on event status.'''
        print(self.ID,self.PreviousEventCheckTime,sep=":")

    def SetTriggerCondition(self,Threshold=0):
        '''Assigns the trigger condition for the satellite. If the event threshold exceeds the assigned trigger condition, a transmission is triggered.
            The trigger condition is intended to be dynamic.'''
        self.TriggerThreshold = Threshold

    def CheckTransmissionConditions(self,TriggerFactor = 1):
        '''Check conditions for successful transmission. These conditions are:
        1. Event threshold check
        2. Radio communication to swarm leader check (TBD)'''
        
        ConditionArray = []
        if TriggerFactor >= self.TriggerThreshold:
            ConditionArray.append(True)
        else:
            ConditionArray.append(False)
        
        if all(ConditionArray) == True:
            self.TriggerTransmission()
        else:
            pass # Do nothing

class Event:
    
    def __init__(self,Name:str=None,Priority:float=0,Threshold:float=0.5,PDF=lambda x:x,RNG:np.random.Generator=None,AdditionalArguments:tuple=None,ThresholdType='FIXED',EventTriggeredFunction=lambda x:x,EventTriggeredArguments:tuple=None):
        '''Initialize the event parameters.
         1. Name -> string: Representative, non-functional
         2. Priority -> float: Event priority
         3. Threshold -> float: Threshold for event trigger probability to be recognized. Limited to [0,1]
         4. PDF -> function: Probability Distribution Function
         5. RNG -> Generator: Seeded RNG function of the form numpy.random.default_rng(#)'''
        self.Name = Name
        self.Priority = Priority
        self.Threshold = max(min(Threshold,1),0)
        self.ThresholdType = ThresholdType
        self.PDF = PDF
        self._RNG = None
        self.EventTriggered = EventTriggeredFunction
        DefaultRNG = np.random.default_rng(58235)
        if RNG is None:
            RNG = DefaultRNG
        self.RNG = RNG # READ-ONLY, assigned once at initialization. Seeded RNG Function in the form numpy.random.default_rng(#) where # is the seed as an integer
        self.LastTriggerTime = 0 # Simulaion time (min) of the previous event trigger
        self.LastCheckTime = 0 # Time at which last event trigger check was performed
        self.AdditionalArguments = AdditionalArguments # Additional arguments for each individual event
        self.EventTriggeredArguments = EventTriggeredArguments # Additional arguments for event triggered function

    @property
    def RNG(self):
        return self._RNG
        
    @RNG.setter
    def RNG(self,RNG):
        if self._RNG == None:
            self._RNG = RNG
        else:
            pass
        return None
    
    def CheckFixedThresholdType(self):
        '''Check if the event has a fixed threshold type.'''
        if self.ThresholdType == 'FIXED':
            return True
        return False

    def CheckEventTrigger(self,Var=0):
        '''Returns True if event triggers with the given parameters, else False.
        The input 'Var' would generally refer to Time but this is not strictly necessary and depends on the underlying PDF.
        How it works:
        If the probability exceeds the threshold then event will trigger, else it will not.'''

        retVal = False
        DeltaVar = Var - self.LastCheckTime
        if isinstance(self.AdditionalArguments,tuple):
            P = self.PDF(Var,DeltaVar,*self.AdditionalArguments)
        else:
            P = self.PDF(Var,DeltaVar)

        #Limit the probablity to the interval [0,1]
        P = max(min(P,1),0)
        if self.CheckFixedThresholdType() == True:
            Threshold = self.Threshold
        else:
            Threshold = self.RNG.random()
            if self.Threshold >= 1:
                Threshold = self.Threshold
        if P >= Threshold:
            retVal = True
            self.LastTriggerTime = Var
            self.State = True
            self.EventTriggerFunction()
        else:
            self.State = False

        self.LastCheckTime = Var

        return retVal
    
    def EventTriggerFunction(self):
        '''Function called when event triggers. This function is used to handle threshold calculation for next trigger if an event is triggered.'''
        if self.EventTriggeredArguments is not None:
            self.Threshold = self.EventTriggered(self.Threshold,*self.EventTriggeredArguments)
        else:
            self.Threshold = self.EventTriggered(self.Threshold)




if __name__ == "__main__":
    pass