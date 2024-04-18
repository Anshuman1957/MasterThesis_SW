from SwarmModel import Swarm
import numpy as np
import scipy.constants

class GroundStation:

    def __init__(self,ID=None):
        '''Initializes the ground station.'''
        self.ID = ID
        self.Swarm = None
        self.Distance = 1.4 # Distance in Astronomical Units
        self.Distance_m = self.Distance * scipy.constants.astronomical_unit # Distance in metres
        self.TransmissionSpeed = scipy.constants.speed_of_light # Speed of transmission
        self.TransmitOneWayTime = self.Distance_m/self.TransmissionSpeed # One-way transmission time
        self.TransmitRoundTripTime = self.TransmitOneWayTime*2 # Round trip transmission time


    def AssignSwarm(self,Swarm:Swarm = None):
        '''Assign the swarm to the ground station. Only 1 swarm can be assigned to a ground station at a time.'''
        if self.Swarm != None:
            print(f"Warning! Re-assigning swarm assigned to Ground Station {self.ID=}")
        self.Swarm = Swarm

    def MoveSwarm(self, Swarm:Swarm, Coordinates:np.array):
        '''Commands the swarm to move from their current position to the given position. The path taken is up to the swarm to decide.'''
        Swarm.MoveSwarm(Coordinates=Coordinates,Relative=False)

    def MoveSwarmRelative(self, Swarm:Swarm, Coordinates:np.array):
        '''Commands the swarm to move from their current position by the given offset. The path taken is up to the swarm to decide.'''
        Swarm.MoveSwarm(Coordinates=Coordinates,Relative=True)

    def GetSwarmHealth(self):
        '''Requests the swarm to return the swarm health status'''
        pass


if __name__ == "__main__":
    pass