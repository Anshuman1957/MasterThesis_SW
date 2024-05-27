from SatelliteModel import Satellite
from EventModel import Event
import numpy as np
import matplotlib.pyplot as plt

class Swarm:

    def __init__(self, SatelliteList:list[Satellite]= [] ):
        '''Initializes the swarm with the given input list of Satellites. The swarm leader is by default set to the first swarm in the provided list.'''

        self.Members = SatelliteList
        self.SwarmSize = len(self.Members)
        self._Leader = []
        self.Leader = SatelliteList[0]
        self.PreviousEventCheckTime = None
        for satellite in self.Members:
            satellite.UpdateSwarmReference(self)

        Satellite_ID_Dict = {}
        for satellite in self.Members:
            Satellite_ID_Dict[satellite.ID] = satellite

        self.Satellite_ID_Dict = Satellite_ID_Dict

        # Additional parameters:
        # Swarm relative distance/position matrices
        # Swarm upkeep time?
        # Swarm housekeeping data (Depends on distance from ground station and mission length)

    def PrintID(self):
        '''Prints the ID of all the satellites in the swarm'''
        for satellite in self.Members:
            print(f'{satellite.ID=}')
        return None

    @property
    def Leader(self):
        '''Getter: Returns the internal leader list'''
        return self._Leader
    
    @Leader.setter
    def Leader(self,Leader):
        '''Adds the input Leader (satellite or list object) to the internal leader list'''
        if isinstance(Leader,Satellite):
            if Leader not in self._Leader:
                self._Leader.append(Leader)
                Leader.Hierarchy = 'Leader'
        elif isinstance(Leader,list):
            for satellite in Leader:
                if satellite not in self._Leader:
                    satellite.Hierarchy = 'Leader'
                    self._Leader.append(satellite)

    @Leader.deleter
    def Leader(self):
        '''Removes all leader satellites from the internal leader list'''
        self._Leader = []

    def AddLeader(self,InSat:Satellite):
        '''Add the given satellite in the swarm to be the leader(s) of the swarm. If the given satellite is not in the swarm, it is added to the swarm.'''
        if InSat not in self.Members:
            self.Members.append(InSat)
        if InSat not in self.Leader:
            self.Leader.append(InSat)
        self.SwarmSize = len(self.Members)
        return None
    
    def UpdateSwarmMembers(self,AddSats:list[Satellite]=None,DelSats:list[Satellite]=None):
        '''Updates (Add or remove) satellites in the swarm followed by updating the related swarm parameters.'''
        localSwarm = self.Members
        # Delete list
        if DelSats is not None:
            tempList = [sat for sat in localSwarm if sat not in DelSats]
            localSwarm = tempList
        # Add list
        if AddSats is not None:
            for Sat in AddSats:
                localSwarm.append(Sat)
        # Update swarm parameters
        self.SwarmSize = len(localSwarm)
        self.Members = localSwarm
        return None
    
    def UpdateSwarmRelativePosition(self,RefLeader:Satellite):
        '''Updates the position of all the swarm satellites with respect to the leader. The leader is treated as the origin'''
        for Satellite in self.Members:
            Satellite.RelativePosition = Satellite.Position - RefLeader.Position # Leader is always at origin with respect to itself
        return None
    
    def PlotCurrentPosition(self, ax=None,fig=None):
        '''Plots the swarm satellites' current absolute positions in a 3D co-ordinate space. The leader is highlighted in red. Raises a warning if any satellite share
        a position'''

        checked = []
        uniquePos = [list(sat.Position) for sat in self.Members if list(sat.Position) not in checked and not checked.append(list(sat.Position))]
        if len(checked) != len(uniquePos):
            raise "Duplicate Positions!"
        else:
            if fig is None:
                fig = plt.figure()
            if ax is None:
                ax = fig.add_subplot(projection='3d')

            for satellite in self.Members:
                ax.scatter(satellite.Position[0],satellite.Position[1],satellite.Position[2],c=('blue' if self.Leader != satellite else 'red'))
            plt.show()

        return (ax,fig)

    def MoveSwarm(self,Coordinates:np.array, Relative=False):
        '''Moves the swarm centered on a swarm satellite to the input location.'''
        pass

    def CheckEvents(self,SimulationTime) -> None:
        '''Update swarm event states based on current simulation time.'''
        
        for satellite in self.Members:
            satellite.UpdateSimulationTime(SimulationTime)
            satelliteEvents = satellite.CheckEventStates()
            satelliteArbitration = satellite.ArbitrateEventStates()
        self.PreviousEventCheckTime = SimulationTime

        return None
    
    def AddCommonEvent(self,Events):
        '''Adds the input Event(s) to all the member satellites, including the Leader. Valid Event input is a list of Event classes or a single Event class.'''
        if isinstance(Events,list):
            EventList = Events
        elif isinstance(Events,Event):
            EventList = [Events]
        else:
            return
        for satellite in self.Members:
            satellite.AddEvents(EventList)

    def UpdateSimulationTime(self,Time:float):
        '''Updates the simulation time for all satellites in the swarm'''
        for satellite in self.Members:
            satellite.UpdateSimulationTime(Time)

    def InitializeInformationMatrix(self):
        '''Trigger initialization of the information matrix of all sub-swarm members'''
        Size = self.SwarmSize
        for satellite in self.Members:
            _ = satellite.InitializeInformationMatrix(Size)
        return None
        
    def SetFollowerTriggerParameter(self,Value:float=1):
        '''Set the epsilon value of each follower satellite to the input value'''
        for satellite in self.Members:
            if satellite not in self.Leader:
                satellite.SetTriggerCondition(Value)
        return None
    
    def SetLeaderTriggerParameter(self,Value:float=1):
        '''Set the epsilon value of each leader satellite to the input value'''
        for satellite in self.Members:
            if satellite in self.Leader:
                satellite.SetTriggerCondition(Value)
        return None

    def SetTriggerParameter(self,Value:float=1):
        '''Set the epsilon value of every satellite to the input value'''
        for satellite in self.Members:
            satellite.SetTriggerCondition(Value)
        return None
    
    def GetSatelliteByID(self,ID:int):
        '''Returns the satellite with the input ID'''
        satellite = self.Satellite_ID_Dict[ID]
        return satellite
    
    def SetMinimumTransmitInterval(self,Interval:int=300):
        '''Sets the minimum transmit interval for each satellite in the swarm'''
        for satellite in self.Members:
            satellite.MinimumTransmitInterval = Interval
        return None
    
    def SetMaximumTransmitInterval(self,Interval:int=60000):
        '''Sets the minimum transmit interval for each satellite in the swarm'''
        for satellite in self.Members:
            satellite.MaximumTransmitInterVal = Interval
        return None
    
    def UpdateIdenticalTransmitFlag(self,State:bool=False):
        '''Updates the flag to enable/disable transmissions containing repeated data of previous transmissions'''
        for satellite in self.Members:
            satellite.IgnoreRepeatData = State
        return None

if __name__ == "__main__":
    x = Satellite(1)
    y = Satellite(2)
    swarm = Swarm([x,y])
    swarm.PrintID()