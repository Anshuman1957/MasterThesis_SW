import scipy.constants
import numpy as np
from SatelliteModel import Satellite
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from OrbitalElements import Orbit

class Asteroid:

    def __init__(self,radius:int,mass:int,period:int=0,position:np.array = np.array([0,0,0]),axis:np.array = np.array([0,0,1])):
        '''Spherical model assumed for now'''
        self.Radius = radius # m
        self.Mass = mass    # kg
        self.RotationPeriod = period # Time taken by asteroid to rotate once in minutes
        self.StandardGravitationalParameter = scipy.constants.gravitational_constant * self.Mass # mu
        self.Gravity = self.StandardGravitationalParameter / (self.Radius**2)
        self.Position = position # Location of the asteroid in 3D coordinates (x,y,z)
        self.Velocity = np.array([0,0,0])
        self.GeoRadius = round(self.ComputeGeostationaryOrbitRadius(),3)
        self.EscapeVelocity = np.sqrt(2*self.Mass * scipy.constants.gravitational_constant / self.Radius)
        #Normalize axis unit vector
        self.Axis = axis / (np.linalg.norm(axis) + 1e-16)

    # What is required here?
    # Asteroid size, heading, velocity, dimensions(?)
    # reflectance(?), can be checked later if time is available
    # Any other info?
        
    def PlotSatelliteOrbit(self,Satellite:Satellite,ax=None,fig=None,PlotOrbit=True,FociPosition=0,PlotInitialPosition=False,debug=False,**kwargs):
        '''Plots the input satellite's orbit around the asteroid considering the satellite's current parameters.
            r = (l^2)/(m^2 * mu) x (1/(1 + e cos (theta)))
            Theta -> Angle that r makes with periapsis of asteroid
            l -> angular momentum of the satellite
            mu -> Standard gravitational parameter
            e -> Eccentricity of the orbit
            m -> Mass of the asteroid
            '''
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(projection='3d')

        Xf = Satellite.Orbit.GetFociValues()[FociPosition][0]
        Yf = Satellite.Orbit.GetFociValues()[FociPosition][1]
        Zf = Satellite.Orbit.GetFociValues()[FociPosition][2]

        X0 = self.Position[0]
        Y0 = self.Position[1]
        Z0 = self.Position[2]

        X0 = X0 - Xf
        Y0 = Y0 - Yf
        Z0 = Z0 - Zf

        self.UpdateSatellitePeriod(Satellite)

        CircleCoordinates = Satellite.Orbit.TranslateOrbit(X0,Y0,Z0)

        Xn = np.array([point[0] for point in CircleCoordinates],dtype=object)
        Yn = np.array([point[1] for point in CircleCoordinates],dtype=object)
        Zn = np.array([point[2] for point in CircleCoordinates],dtype=object)

        #OrbitCoordinates = np.array(list(zip(Xn,Yn,Zn)))
        if PlotOrbit == True:
            if 'SatelliteOrbitColor' in kwargs.keys():
                color = kwargs['SatelliteOrbitColor']
            else:
                color = 'b'
            if 'OrbitAlpha' in kwargs.keys():
                alpha = kwargs['OrbitAlpha']
            else:
                alpha = 0.2
            ax.plot(Xn,Yn,Zn,alpha=alpha,color=color)
        if PlotInitialPosition == True:
            if 'SatelliteMarkerColor' in kwargs.keys():
                color = kwargs['SatelliteMarkerColor']
            else:
                color = 'b'
            ax.scatter(Xn[0],Yn[0],Zn[0],c=color)
        ax.set_aspect('equal',anchor='C')
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_zlabel('Z-axis')
        #ax.set_xlim([-30000, 30000])
        #ax.set_ylim([-30000, 30000])
        #ax.set_zlim([-30000, 30000])

        return (ax,fig)

    def ComputeSatelliteAltitude(self,OrbitalPeriod:int,ReturnAltitude=True):
        '''Computes the value of the semi-major axis or altitude in kilometres given the required orbital period using the formula:
        T = 2*pi * sqrt(a^3/mu)
        Where T -> Orbital Period in minutes
        a -> semi-major axis in kilometres
        mu -> Standard Gravitational Parameter of the asteroid
        '''

        T = OrbitalPeriod
        mu = self.StandardGravitationalParameter
        a = np.cbrt(mu*(T**2)/(4*(np.pi**2)))

        if ReturnAltitude == True:
            return round((a-self.Radius),2)
        else:
            return round(a,2)
    
    def ComputeSatellitePeriod(self,OrbitAltitude:int,isAltitude=False):
        '''Computes the period of rotation for the input orbital altitude or radius. The formula used is:
        T = 2*pi * sqrt(a^3/mu)
        Where T -> Orbital Period in minutes
        a -> semi-major axis in kilometres
        mu -> Standard Gravitational Parameter of the asteroid
        '''
        if isAltitude == True:
            a = OrbitAltitude + self.Radius
        else:
            a = OrbitAltitude
        mu = self.StandardGravitationalParameter
        T = 2*np.pi*np.sqrt((a**3)/mu)
        T = round(T,2)

        return T
    
    def UpdateSatellitePeriod(self,Satellite:Satellite):
        '''Update the input satellite's period with it's currently set radius'''
        
        SatelliteOrbitPositions = Satellite.Orbit.GetOrbitValues()
        Satellite.Orbit.OrbitalPeriod = self.ComputeSatellitePeriod(Satellite.Orbit.Radius)
        Satellite.Orbit.StoreOrbitValues(SatelliteOrbitPositions)
        return None

    def PlotAsteroid(self,ax=None,fig=None,PlotAxis=False):
        '''Plots a spherical model of the asteroid and returns the axis and figure objects'''
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(projection='3d')

        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = ((self.Radius/1) * (np.cos(u)*np.sin(v))) + self.Position[0]
        y = ((self.Radius/1) * (np.sin(u)*np.sin(v))) + self.Position[1]
        z = ((self.Radius/1) * (np.cos(v))) + self.Position[2]
        ax.plot_surface(x, y, z, linewidth=0, color='brown', alpha = 0.5)
        x = self.Position[0]
        y = self.Position[1]
        z = self.Position[2]
        u = self.Axis[0]
        v = self.Axis[1]
        w = self.Axis[2]
        if PlotAxis == True:
            ax.quiver(x,y,z,u,v,w,pivot='middle',normalize=False,length = self.Radius*2.4)
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_zlabel('Z-axis')
        return (ax,fig)
    
    def ComputeGeostationaryOrbitRadius(self):
        '''Computes the geostationary orbital radius using the below formula:
        r = Cuberoot(G * M * T^2)/(4 * pi^2)
        Where:
        r -> Orbit radius
        G -> Gravitational constant
        M -> Mass of the asteroid
        T -> Rotation period of the asteroid'''

        return np.cbrt((scipy.constants.gravitational_constant * self.Mass * (self.RotationPeriod*60)**2)/(4*(np.pi**2))) 

class AsteroidField:

    def __init__(self,AsteroidList:list[Asteroid]):
        '''Collection of asteroids'''
        self.Asteroids = AsteroidList
    
    def PlotAsteroidField(self,ax=None,fig=None,PlotAxis=True):
        '''Plots the entire asteroid field in 3-D coordinates'''

        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111,projection='3d')
        

        def AsteroidData(asteroid:Asteroid):
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x = ((asteroid.Radius/1000) * (np.cos(u)*np.sin(v))) + asteroid.Position[0]
            y = ((asteroid.Radius/1000) * (np.sin(u)*np.sin(v))) + asteroid.Position[1]
            z = ((asteroid.Radius/1000) * (np.cos(v))) + asteroid.Position[2]
            return x,y,z
        
        def PlotAsteroidData(data,ax=None,fig=None,PlotAxis=True):
            if fig is None:
                fig = plt.figure(figsize=(6,6))
            if ax is None:
                ax = fig.add_subplot(111,projection='3d')
            for index,ast in enumerate(data):
                x,y,z = ast[0],ast[1],ast[2]
                ax.plot_surface(x,y,z,linewidth=0, color='brown',alpha=0.5)
                if PlotAxis == True:
                    for asteroid in self.Asteroids:
                        x = asteroid.Position[0]
                        y = asteroid.Position[1]
                        z = asteroid.Position[2]
                        u = asteroid.Axis[0]
                        v = asteroid.Axis[1]
                        w = asteroid.Axis[2]
                        ax.quiver(x,y,z,u,v,w,pivot='middle',normalize=False,length = np.sqrt(asteroid.Radius//20))
            ax.set_aspect('equal')
            ax.set_xlabel('X-axis')
            ax.set_ylabel('Y-axis')
            ax.set_zlabel('Z-axis')
            plt.show()

        data = [AsteroidData(asteroid) for asteroid in self.Asteroids]
        PlotAsteroidData(data,ax,fig,PlotAxis=PlotAxis)


        
        return (ax,fig,data)


Eros = Asteroid(8420,6.687*(10**15)) # 433 Eros

if __name__ =="__main__":
    print(f"{Eros.Gravity=}")