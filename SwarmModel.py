from SatelliteModel import Satellite, Event
import numpy as np
import matplotlib.pyplot as plt
from OrbitalElements import Orbit,Rotation

class Swarm:

    def __init__(self, SatelliteList:list[Satellite]= [] ):
        '''Initializes the swarm with the given input list of Satellites. The swarm leader is by default set to the first swarm in the provided list.'''

        self.Members = SatelliteList
        self.SwarmSize = len(self.Members)
        self.Leader = SatelliteList[0]
        self.PreviousEventCheckTime = None

        # Additional parameters:
        # Swarm size (self.SwarmSize)
        # Swarm relative distance/position matrices
        # Swarm upkeep time?
        # Swarm housekeeping data (Depends on distance from ground station and mission length)
        # Packet losses (Noise in transmission)

    def PrintID(self):
        '''Prints the ID of all the satellites in the swarm'''
        for satellite in self.Members:
            print(f'{satellite.ID=}')
        return None


    def AssignLeader(self,InSat:Satellite):
        '''Assigns the given satellite in the swarm to be the leader of the swarm. If the given satellite is not in the swarm, it is added to the swarm.'''
        if InSat not in self.Members:
            self.Members.append(InSat)
        self.Leader = InSat
        return None
    
    def UpdateSwarmMembers(self,AddSats:list[Satellite]=None,DelSats:list[Satellite]=None):
        '''Updates (Add or remove) satellites in the swarm followed by updating the related swarm parameters. Resets the leader to the 0th Satellite if the leader was removed.'''
        localSwarm = self.Members
        PrevLeader = self.Leader
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
        if PrevLeader not in self.Members:
            self.Leader = self.Members[0]
        return None
    
    def UpdateSwarmRelativePosition(self):
        '''Updates the position of all the swarm satellites with respect to the leader. The leader is treated as the origin'''
        for Satellite in self.Members:
            Satellite.RelativePosition = Satellite.Position - self.Leader.Position # Leader is always at origin with respect to itself
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

    def CheckEvents(self,SimulationTime):
        '''Checks swarm events based on current simulation time.'''
        
        for satellite in self.Members:
            satelliteEvents = satellite.CheckEventStates(SimulationTime)
            satelliteArbitration = satellite.ArbitrateEventStates()
        self.PreviousEventCheckTime = SimulationTime
    
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
        

if __name__ == "__main__":
    x = Satellite(1)
    y = Satellite(2)
    swarm = Swarm([x,y])
    swarm.PrintID()