import numpy as np
import matplotlib.pyplot as plt
from OrbitalElements import Orbit
import FunctionLibrary
import copy
from EventModel import Event

class Satellite:


    def __init__(self,ID:int = 0,InitialPosition:np.array = np.array([0,0,0]),InitialRelativePosition:np.array = np.array([0,0,0])):
        '''Initializes the Satellite with the following information:
        1. Identification Number for the satellite (Default 0)
        2. Co-ordinate position (Defaults to 3D co-ordinate geometry. Can choose quaternion representation)
        3. Initializes the parameter values (TBD)'''
        self.ID = ID # Identifier number for the satellite in the swarm
        self.Mass = 1 # kg, Pico-satellite
        self.InitializePosition(InitialPosition) # Initialize the position of the satellite
        self.InitializeRelativePosition(InitialRelativePosition) # Initialize the relative position of the satellite
        self.Velocity = np.array([0,0,0]) # Velocity of the satellite in (x,y,z) in m/s
        self.Events = [] # Initialize an empty array of events attached to this satellite
        self.EventStates = [] # Current state of events (boolean)
        self.PreviousEventCheckTime = None # Simulation time in minutes when the last event check occured
        self._TriggerThreshold = 0
        if self.ID != 0: # Followers
            self.SetTriggerCondition(0.26) # 1/sqrt(SwarmSize)
            #self.SetTriggerCondition(0.577)
        else:
            self.SetTriggerCondition(1)
        self.SimulationTime = 0
        self.ClearPositionHistory()
        self.ClearRelativePositionHistory()
        self.InformationMatrixHistory = []
        self.UpdateTransmissionParameters()
        self.ReceivedTransmissions = {}
        self.Hierarchy = 'Follower'
        self.MIN_DISTANCE_FACTOR = 4.605
        self.LastTransmissionData = {'Event Information':[]}
        self.TransmissionHistory = []

        # EXPERIMENTAL
        self.LastTransmitTime = -1000
        self.MinimumTransmitInterval = 300
        self.MaximumTransmitInterVal = 60000
        self.TX_THIS_TIME = False # Used for plotting
        self.U_DIFF_History = []
        self.IgnoreRepeatData = False
        #self.Type = # Satellite role in the heterogeneous swarm

    @property
    def LastTransmitTime(self):
        return self._LastTransmitTime
    
    @LastTransmitTime.setter
    def LastTransmitTime(self,Time:int):
        self._LastTransmitTime = Time

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
    
    def UpdateSwarmReference(self,InputSwarm:'Swarm'=None): # type: ignore
        '''Circular reference to the parent swarm relevant for transmission trigger'''
        self.SwarmReference = InputSwarm

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

    def UpdateSimulationTime(self,Time:float):
        '''Updates the stored simulation time. To be called during the time super-loop in main'''
        self.PreviousEventCheckTime = self.SimulationTime
        self.SimulationTime = Time
    
    def AddEvents(self,events):
        '''Adds the input Event(s) to the satellite event handler. Sets the initial event state to False.'''
        if isinstance(events,list):
            for event in events:
                if event not in self.Events:
                    self.Events.append(event)
                    self.EventStates.append(event.State)
                    self.LastTransmissionData['Event Information'].append(event.State)
        elif isinstance(events,Event):
            if events not in self.Events:
                self.Events.append(events)
                self.EventStates.append(events.State)
                self.LastTransmissionData['Event Information'].append(events.State)
        else:
            raise TypeError
        
        
        return None

    def ClearEvents(self):
        '''Clears all attached events to the satellite'''
        self.Events = []
        self.EventStates = []

        return None

    def CheckEventStates(self):
        '''Check if any events trigger based on simulation time. Returns a list of boolean values corresponding to the event trigger state.'''
        T = self.SimulationTime
        self.PreviousEventCheckTime = T
        retVal = []
        for index,event in enumerate(self.Events):
            EventState = event.CheckEventTrigger(T)
            retVal.append(EventState)
            self.EventStates[index] = EventState
        self.CurrentEventStates = retVal
        return retVal
    
    def ArbitrateEventStates(self,Events:list[Event]=None,EventStates=None):
        '''Compute any external action required based on the event states. Call TriggerTransmission() if event(s) exceed defined action threshold.'''
        if Events is None:
            Events = self.Events
        if EventStates is None:
            EventStates = self.CurrentEventStates

        T = self.SimulationTime
        Priorities = np.array([event.Priority(T) for event in Events])
        EventStatesInt = np.array(EventStates,dtype=bool)
        EventStatesInverted = np.invert(EventStatesInt)
        X_Matrix = Priorities*EventStatesInverted
        TriggerFactor = np.linalg.norm(X_Matrix,ord=2)#sum(Priorities) - np.linalg.norm(X_Matrix,ord=2)
        TriggerFactor = round(TriggerFactor,4)
        if self.CheckTransmissionConditions(TriggerFactor):
            #self.SetTriggerCondition(Threshold=1.00,Type='MULTIPLY')
            pass

        return None

    
    def TriggerTransmission(self,ForceTX=False) -> None:
        '''Trigger a (short-range radio) transmission broadcast intended for the leader of the swarm based on event status.'''
        #Construct the transmission information dictionary
        TransmitInfo = {}
        TransmitInfo['ID'] = self.ID
        TransmitInfo['Time'] = self.SimulationTime
        TransmitInfo['Position'] = self.Position
        EventStates = list(np.array(self.EventStates,copy=True))
        TransmitInfo['Event Information'] = EventStates # Update with specific event information
        TransmitInfo['Events'] = [event.Name for event in self.Events]
        TransmitInfo['Event Decay Time'] = [event.HoldTime for event in self.Events]
        LastTXEventStates = self.LastTransmissionData['Event Information']
        
        if (EventStates == LastTXEventStates) and (ForceTX == False) and (self.IgnoreRepeatData == True): # If event info data is same as last time and it is not a forced Tx, abort transmission
            return None
        self.TX_THIS_TIME = True # Used for plotting purposes only
        self.TransmissionHistory.append(TransmitInfo)
        self.LastTransmitTime = copy.deepcopy(self.SimulationTime)
        #if ForceTX == False: # Do not update Ucap if it is a forced transmission
        self.U_CAP = np.array(self.InformationMatrix,copy=True)

        for satellite in self.SwarmReference.Members:
            _ = satellite.ReceiveTransmission(TransmitInfo)


    def ReceiveTransmission(self,TXData:dict=None) -> None:
        '''Receive a transmission broadcast'''
        
        Distance = FunctionLibrary.P2Length(TXData['Position'],self.Position)
        
        if not isinstance(TXData,dict): # Ignore incorrect TXData format
            pass
        
        elif TXData['ID'] == self.ID: # Ignore receiving your own transmission
            self.LastTransmissionData = TXData
        
        elif Distance > self.ReceiverRange: # Ignore transmission if out of reception range
            pass
        
        else:
        # Transmission received
            self.ReceivedTransmissions[TXData['ID']] = TXData
        

        # Consideration to use the received event data to compute information factor for an amount of time

        return None

        

    def SetTriggerCondition(self,Threshold=0,Type=None):
        '''Assigns the trigger condition for the satellite. If the event threshold exceeds the assigned trigger condition, a transmission is triggered.
            The trigger condition is intended to be dynamic.'''
        if Type == 'ADD':
            self.TriggerThreshold = self.TriggerThreshold + Threshold
        elif Type == 'MULTIPLY':
            self.TriggerThreshold = self.TriggerThreshold * Threshold
        elif Type == None:
            self.TriggerThreshold = Threshold
    

    def CheckTransmissionConditions(self,TriggerFactor = 1):
        '''Check conditions for successful transmission. These conditions are:
        1. Event threshold check
        2. Radio communication to swarm leader check (TBD)
        3. Minimum time since last transmission check
        4. Maximum time since last transmission check
        5. Change in data. Check is performed in TriggerTransmission() function as the message is constructed there.'''

        self.TX_THIS_TIME = False
        retVal = None
        ConditionArray = []

        # Min/Max time since last transmission check
        TimeSinceLastTransmit = self.SimulationTime - self.LastTransmitTime
        if (TimeSinceLastTransmit > self.MinimumTransmitInterval):
            ConditionArray.append(True)
        else:
            ConditionArray.append(False)

        # Check || U - u ||_2 <= Epsilon*TriggerFactor, Norm-2 of current and previous information matrix
        InfoMatrixLocal = np.array(self.InformationMatrix,copy=True)
        U_DIFF = np.linalg.norm(self.U_CAP - InfoMatrixLocal,ord=2)
        self.U_DIFF_History.append(U_DIFF)
        U_DIFF = round(U_DIFF,4)
        TriggerThreshold = round(TriggerFactor*self.TriggerThreshold,4)
        if U_DIFF >= (TriggerThreshold):
            ConditionArray.append(True)
        else:
            ConditionArray.append(False)
        if (TimeSinceLastTransmit >= self.MaximumTransmitInterVal):
            self.TriggerTransmission(ForceTX=True)
            retVal = True
        elif (all(ConditionArray) == True):
            self.TriggerTransmission()
            retVal = True
        else:
            pass # Do nothing

        return retVal
    def ComputeOrbitParameters(self, OrbitPeriod:int,OrbitFoci:int=0) -> None:
        '''Computes the orbital time parameters for the input asteroid object.'''
        OrbitPositions = self.Orbit.GetOrbitValues()
        if OrbitFoci > 1:
            OrbitFoci = 0
        Focus = self.Orbit.GetFociValues()[OrbitFoci]
        T = OrbitPeriod
        n = self.Orbit.OrbitalPoints
        OrbitFraction = []
        a = self.Orbit.Radius
        b = a * np.sqrt(1-(self.Orbit.Eccentricity**2))
        totalArea = round(np.pi * a * b,2)
        intA = 0
        for index,point in enumerate(OrbitPositions):
            if index == 0:
                point1 = point
                point2 = point
            else:
                point1 = point
                point2 = OrbitPositions[(index-1)%n]
            # Compute area of triangle
            dA = FunctionLibrary.P3Area(point1,point2,Focus)
            # Compute fraction of area covered so far
            intA = intA + dA
            # Convert area to time
            intT = round(intA / totalArea,4)
            # Convert time to orbit ratio
            dT = round(intT * T,4)
            # Append to OrbitFraction array
            OrbitFraction.append(dT)

        OrbitFraction.append(T) # The first and last time values are for the starting position of the orbit
        self.Orbit.StoreTimeValues(OrbitFraction,T)

        return None

    def InitializeInformationMatrix(self,Size:int) -> None:
        '''Initialize the Information matrix as an n-by-n matrix filled with zeros'''
        self.InformationMatrix = np.zeros((Size,Size))
        self.InformationMatrixPrev = self.InformationMatrix
        InformationMatrix = np.array(self.InformationMatrix,copy=True)
        self.U_CAP = np.linalg.norm(InformationMatrix,ord=2)
        self.InformationMatrixHistory = []
        self.TransmissionHistory = []
        del InformationMatrix

        return None
    
    def StoreInformationMatrix(self) -> None:
        '''Updates the previous Information Matrix value with the current value'''
        InformationMatrix = np.array(self.InformationMatrix,copy=True)
        self.InformationMatrixHistory.append(InformationMatrix)
        self.InformationMatrixPrev = InformationMatrix
        del InformationMatrix

        return None
    
    def CheckInformationMatrixJacobian(self,dt:int=1) -> None:
        '''Computes the Jacobian of the information matrix by taking the difference of the matrix with the previous iteration matrix'''
        Jacobian = (self.InformationMatrix - self.InformationMatrixPrev)/dt
        # Check || U - u ||2 <= Epsilon, Norm-2 of current and previous information matrix
        uCap = np.array(self.InformationMatrix)
        u = np.array(self.InformationMatrixPrev)
        Diff = np.linalg.norm(uCap - u,ord=2)
        #self.U_CAP = Diff

        return None
    
    def UpdateSelfInformationFactor(self,Alpha:float=0.5,Beta:float=0.5,Time:float=0,ExpectedPosition:np.array=None,Noise:np.array=None,dMax:float = 20):
        '''Computes,stores and returns the self-information factor of the satellite. If the ID of the satellite is 'x', then this self-information factor is \
            a_(xx)'''
        MINIMUM_INFORMATION_FACTOR = self.MIN_DISTANCE_FACTOR # Derived by solving the following for x: e^x = 1%
        CurrentPosition = self.Position
        if ExpectedPosition is not None:
            ExpectedPosition = ExpectedPosition # Does nothing tbh, as long as it is a numpy array representing a coordinate
        else:
            if Noise is None:
                ExpectedPosition = CurrentPosition
            else:
                ExpectedPosition = CurrentPosition + Noise
        FactorA = FunctionLibrary.InfoFactorA(Alpha=Alpha,Time=Time,EventList=self.Events,EventStates=self.EventStates,Estimate=False)
        FactorB = FunctionLibrary.InfoFactorB(Beta,CurrentPosition,ExpectedPosition,dMax,MINIMUM_INFORMATION_FACTOR)
        InformationFactor = round(FactorA + FactorB,4)
        self.InformationFactor = InformationFactor
        self.InformationMatrix[self.ID][self.ID] = InformationFactor

        return InformationFactor
    
    def UpdateSelfCrossInformationFactor(self,InSat:'Satellite',Time:float=0,Alpha=0.5,Beta=0.5,dMax = 20,PosNoise = None,EventsKnown = False,Estimate=True, PositionKnown = False,delTMax:int = 0,Invert=False):
        '''Computes, stores, and returns the cross-information factor of the satellite. If the ID of the (self) satellite is 'x' and cross satellite is 'y', then this information factor is a_(xy) which indicates what satellite 'x' knows about satellite 'y'. Also returns the self and input satellite's IDs.'''
        # Estimate or Measure based on last received transmission and event decay time:
        EstimateLoc = Estimate
        Measure = None
        # Check if transmission received at all
        if InSat.ID not in self.ReceivedTransmissions.keys():
            pass
        else:
            # Construct decay time matrix
            currentTime = self.SimulationTime
            lastRXTime = self.ReceivedTransmissions[InSat.ID]['Time']
            TimeSinceLastRx = currentTime - lastRXTime
            DecayArray = np.array(self.ReceivedTransmissions[InSat.ID]['Event Decay Time'])
            StateArray = np.array(self.ReceivedTransmissions[InSat.ID]['Event Information'])
            UpdatedEventArray = DecayArray * StateArray # This array contains the hold time for each event
            EstimateLoc = TimeSinceLastRx >= UpdatedEventArray
            if any(EstimateLoc):
                Estimate = False
                Measure = False
        # Alpha term calculation
        if Invert == True:
            InSatRef = InSat
            InSat = self
        if EventsKnown == False:
            term1 = 0
        else:
            term1 = FunctionLibrary.InfoFactorA(Alpha=Alpha,Time=Time,EventList=InSat.Events,Estimate=Estimate,Measure=Measure)

        # Beta term calculation
        if PositionKnown == False:
            term2 = 0
        else:
            InSatPos = InSat.Position
            if PosNoise is None:
                InSatExpectedPosition = InSatPos
            else:
                InSatExpectedPosition = InSatPos + PosNoise
            term2 = FunctionLibrary.InfoFactorB(Beta = Beta,P1=InSat.Position,P2=InSatExpectedPosition,P_MAX=dMax,MIN_FACTOR=InSat.MIN_DISTANCE_FACTOR)
                
        CrossFactor = round(term1 + term2,4)
        if Invert == False:
            self.InformationMatrix[InSat.ID][self.ID] = CrossFactor
        else:
            self.InformationMatrix[self.ID][InSatRef.ID] = CrossFactor

        return (CrossFactor,self.ID,InSat.ID)
    
    def UpdateOtherCrossInformationFactor(self,InSat:'Satellite',InSat2:'Satellite',Time:float=0,Alpha=0.5,Beta=0.5,dMax=20,PosNoise=None,EventsKnown=False,Estimate=True,PositionKnown=False,delTMax:int=0):
        '''Computes, stores, and returns the cross-information factor of the satellite. If the ID of the (self) satellite is 'x' and cross satellite is 'y', then this information factor is a_(yx) which indicates what satellite 'y' knows about satellite 'x'. Also returns the self and input satellite's IDs.'''
        # Alpha term calculation
        if EventsKnown == False:
            term1 = 0
        else:
            term1 = FunctionLibrary.InfoFactorA(Alpha=Alpha,Time=Time,EventList=InSat2.Events,Estimate=Estimate)
        
        # Beta term calculation
        if PositionKnown == False:
            term2 = 0
        else:
            SatPos = InSat2.Position
            if PosNoise is None:
                SatExpectedPosition = SatPos
            else:
                SatExpectedPosition = SatPos + PosNoise
            term2 = FunctionLibrary.InfoFactorB(Beta=Beta,P1=InSat2.Position,P2=SatExpectedPosition,P_MAX=dMax,MIN_FACTOR=self.MIN_DISTANCE_FACTOR)

        CrossFactor = term1 + term2
        self.InformationMatrix[InSat.ID][InSat2.ID] = CrossFactor
        return (CrossFactor,InSat.ID,InSat2.ID)

    def UpdateOtherSelfInformationFactor(self,InSat:'Satellite',Beta:float=1,Noise:np.array=None,NoiseFactor:np.array=None,dMax:float=20,Estimate:bool=True):
        '''Estimate the self-information factor of the other input satellite.
        Based on the last transmission received from the target satellite, the target self-information factor is scaled negatively based on the input Noise Factor arg.'''
        CurrentPosition = InSat.Position
        if InSat.ID in self.ReceivedTransmissions.keys():
            TimeSinceLastReception = self.SimulationTime - self.ReceivedTransmissions[InSat.ID]['Time']
        else:
            TimeSinceLastReception = self.SimulationTime
        DecayFactor = 1000
        NoiseFactor = NoiseFactor * (TimeSinceLastReception/DecayFactor)
        # self.ReceivedTransmissions
        if Estimate == True:
            ExpectedPosition = CurrentPosition + (Noise * NoiseFactor)
        else:
            ExpectedPosition = CurrentPosition + Noise
        IF_A = FunctionLibrary.InfoFactorA(Alpha=Beta,Time=0,EventList=InSat.Events,EventStates=InSat.EventStates,Estimate=False,Measure=True)
        IF_B = FunctionLibrary.InfoFactorB(Beta=Beta,P1=CurrentPosition,P2=ExpectedPosition,P_MAX=dMax,MIN_FACTOR=InSat.MIN_DISTANCE_FACTOR)
        InformationFactor = IF_A + IF_B
        self.InformationMatrix[InSat.ID][InSat.ID] = InformationFactor

        return (InformationFactor,InSat.ID)


    def UpdateTransmissionParameters(self,Power=10,Gain=20,FoV=360,Range=5000):
        '''Updates the radio transmission parameters'''
        self.TransmitterPower = Power # Watts
        self.TransmitterGain = Gain # dB
        self.TransmitterFoV = FoV # degrees
        self.TransmitterRange = Range # km
        self.ReceiverRange = Range # km


if __name__ == "__main__":
    pass