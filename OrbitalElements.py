import scipy.constants
import numpy as np
#from SatelliteModel import Satellite
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from collections import deque


def GenerateCircleCoordinates(Radius,Eccentricity,NumPoints):
    a = Radius
    b = a* np.sqrt(1-(Eccentricity**2))
    return np.array([(round(np.cos(2*np.pi/NumPoints*x)*a,4),round(np.sin(2*np.pi/NumPoints*x)*b,4),0) for x in range(0,NumPoints)])

class Rotation:
    
    def __init__(self,Axis:str='x',Angle:int=0):
        '''Initializes the rotation object'''

        self.Axis = Axis #x,y or z axis
        self.Angle = Angle # Angle in degrees
        self.RotationMatrix = self.ComputeRotationMatrix()

        @property
        def Axis(self):
            return self._Axis
        
        @Axis.setter
        def Axis(self,Axis='x'):
            if Axis in ['x','y','z']:
                self._Axis = Axis
            else:
                self._Axis = 'x'
    
    def __str__(self):
        return f"<{self.Axis} -> {self.Angle} degrees>"
    
    def ComputeRotationMatrix(self):
        '''Computes the rotation matrix for the given axis/angle combination and stores it as a 3x3 numpy array'''
        t = np.deg2rad(self.Angle) # Convert degrees to radians because np.sin and np.cos default argument is in radians
        sin = np.sin
        cos = np.cos
        if self.Axis == 'x':
            Matrix = np.array([[1,0,0],
                               [0,cos(t),-sin(t)],
                               [0,sin(t),cos(t)]],dtype=object)
        elif self.Axis == 'y':
            Matrix = np.array([[cos(t),0,sin(t)],
                                [0,1,0],
                                [-sin(t),0,cos(t)]],dtype=object)
        elif self.Axis == 'z':
            Matrix = np.array([[cos(t),-sin(t),0],
                               [sin(t),cos(t),0],
                               [0,0,1]],dtype=object)
        return Matrix
    
class Orbit:

    def __init__(self,radius:int = 1, e:int=0,Rotations:list[Rotation]=[Rotation('x',0),Rotation('y',0)]):
        '''Initializes the orbit parameters. The orbit is always considered to be centered at the origin (0,0,0). Handling the orbit center should be done by
            the plotting functions by adding the central mass' co-ordinates to the orbit co-ordinates. Alternatively, use the TranslateOrbit() class function to
            generate a set of translated orbits.'''

        # Only circular and elliptical orbits are considered at this time. Hyperbolic (e>1) and parabolic (e=1) orbits are not considered.
        if e < 0:
            self.Eccentricity = 0
        elif e >= 1:
            self.Eccentricity = 0.99
        else:
            self.Eccentricity = e

        self.Radius = radius # Radius in kilometers
        #self.OrbitalPeriod = T # Orbital period expressed in minutes
        self.OrbitalPeriod = 1 # In minutes, to be updated when defining for a specific body
        self.BasePlane = 'xy' # Initial plane is fixed to X-Y plane. It is not meant to be changed.
        self.Rotations = Rotations
        self.OrbitalPoints = 360 # Number of positions to store in the orbit
        
        OrbitCoordinates = GenerateCircleCoordinates(self.Radius,self.Eccentricity,self.OrbitalPoints)
        _ = self.StoreOrbitValues(OrbitCoordinates)
        _ = self.RotateOrbit()

        self.OrbitPeriodRemainder = 0 # Variable to hold the remainder of division when plotting orbital position with time

    @property
    def BasePlane(self):
        '''Property BasePlane is not used at this time'''
        return self._BasePlane
    
    @BasePlane.setter
    def BasePlane(self,NewPlane='xy'):
        '''Sets the base plane for the orbit in the xy,yz or xz axis. If input is invalid, defaults to 'xy' plane'''
        if NewPlane in ['xy','yx']:
            self._BasePlane = 'xy'
        elif NewPlane in ['yz','zy']:
            self._BasePlane = 'yz'
        elif NewPlane in ['xz','zx']:
            self._BasePlane = 'xz'
        else:
            self._BasePlane = 'xy'

    @property
    def Rotations(self):
        return self._Rotations
    
    @Rotations.setter
    def Rotations(self,Rotations:list[Rotation]):
        self._Rotations = Rotations

    @property
    def OrbitalPeriod(self):
        return self._OrbitalPeriod
    
    @OrbitalPeriod.setter
    def OrbitalPeriod(self,Period):
        self._OrbitalPeriod = Period # In minutes
    

    def StoreOrbitValues(self,OrbitalPositions:np.array):
        '''Stores the orbital positions within the Orbit object for later use (mainly plotting)'''
        self._OrbitPoints = OrbitalPositions
        self.PeriodResolution = round(self.OrbitalPeriod / self.OrbitalPoints,2)
        return None
    
    def GetOrbitValues(self):
        '''Returns the stored orbital position values'''
        return self._OrbitPoints

    def ShiftOrbit(self,n=1,reverse=False):
        '''Rotate the orbital points by the input n value. Used to set the initial position in the orbit. 'reverse' input is currently unused as retrograde orbits are not considered.'''
        OrbitCoordinates = self.GetOrbitValues()
        ShiftedCoordinates = np.roll(OrbitCoordinates,3*n) # n is multiplied by 3 because Coordinates is stored as a nx3 matrix (3 dimensions)
        self.StoreOrbitValues(ShiftedCoordinates)

    
    def RotateOrbit(self):
        '''Rotates the orbit around the origin based on the attached Rotation objects and stores the values.
           A maximum of 2 rotations will be processed in order of their existence in self.Rotations'''

        rotations = self.Rotations # Private member, assigned to a local variable
        if len(rotations) == 2: 
            CombinedRotationMatrix = np.dot(self.Rotations[0].RotationMatrix,self.Rotations[1].RotationMatrix)
        else:
            CombinedRotationMatrix = self.Rotations[0].RotationMatrix

        OrbitCoordinates = self.GetOrbitValues()
        RotatedCoordinates = np.dot(OrbitCoordinates,CombinedRotationMatrix)
        _ = self.StoreOrbitValues(RotatedCoordinates)

        return None
    
    def ResetRotations(self):
        '''Resets the orbit coordinates to if no rotations were performed, i.e, on the xy plane.'''
        OrbitCoordinates = GenerateCircleCoordinates(self.Radius,self.Eccentricity,self.OrbitalPoints)
        _ = self.StoreOrbitValues(OrbitCoordinates)


    def TranslateOrbit(self,X0=0,Y0=0,Z0=0):
        '''Returns the orbit coordinates centered at the input X0,Y0,Z0 coordinates. This does not modify the orbit's coordinates centered around the origin.'''
        Coordinates = self.GetOrbitValues()
        X = np.array([point[0] for point in Coordinates])
        Y = np.array([point[1] for point in Coordinates])
        Z = np.array([point[2] for point in Coordinates])
        Xn = X + X0
        Yn = Y + Y0
        Zn = Z + Z0
        CoordinatesTranslated = np.array(list(zip(Xn,Yn,Zn)))
        # There is probably a numpy function to do this in-place faster which I could try to explore, later, possibly

        return CoordinatesTranslated

        