import numpy as np

class Event:
    
    def __init__(self,Name:str=None,Priority:float=lambda x:1,Threshold:float=lambda x:0.5,PDF=lambda x:x,RNG:np.random.Generator=None,AdditionalArguments:tuple=None,ThresholdType='FIXED',InitialState:bool=False, HoldTime:float = -1,SingleTrigger:bool=False):
        '''Initialize the event parameters.
         1. Name -> string: Representative, non-functional
         2. Priority -> float: Event priority
         3. Threshold -> float: Threshold for event trigger probability to be recognized. Limited to [0,1]
         4. PDF -> function: Probability Distribution Function
         5. RNG -> Generator: Seeded RNG function of the form numpy.random.default_rng(#)
         6. AdditionalArguments -> Tuple: A set of additional arguments (optionally) used by the PDF function
         7. ThresholdType -> string: 'FIXED' and other.
         8. InitialState -> bool: Initial boolean state of the event. Default False
         9. HoldTime -> int: Hold time for the event state when it is triggered. Negative indicates until reset.'''
        self.Name = Name
        self.Priority = Priority
        self.PriorityOriginal = Priority
        #self.Threshold = max(min(Threshold,1),0)
        self.Threshold = Threshold
        self.ThresholdType = ThresholdType
        self.PDF = PDF
        self._RNG = None
        self.SingleTrigger = SingleTrigger
        self.State = InitialState
        self.TriggerFlag = False
        DefaultRNG = np.random.default_rng(58235)
        if RNG is None:
            RNG = DefaultRNG
        self.RNG = RNG # READ-ONLY, assigned once at initialization. Seeded RNG Function in the form numpy.random.default_rng(#) where # is the seed as an integer
        self.LastTriggerTime = 0 # Simulaion time (min) of the previous event trigger
        self.LastCheckTime = 0 # Time at which last event trigger check was performed
        self.AdditionalArguments = AdditionalArguments # Additional arguments for each individual event
        self.Error = False # Value to store any event fault which prevents the event from being recognized or measured
        self.HoldTime = HoldTime # minutes, Denotes how long an Event 'TRUE' state should be held when the event is triggered. A negative value indicates indefinitely.

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
    
    def ResetEventState(self) -> None:
        '''Resets the event state. Required when the hold time is indefinite.'''
        self.State = False

        return None

    def CheckEventTrigger(self,Var:float=0,LastCheckTime=None,Estimate:bool=False,Measure:bool=False):
        '''Returns True if event triggers with the given parameters, else False. Can also return 0.5 in the event of an error causing InfoFactorA = 0 for the event.
        The input 'Var' would generally refer to Time but this is not strictly necessary and depends on the underlying PDF.
        How it works:
        If the probability exceeds the threshold then event will trigger, else it will not.'''
        
        if LastCheckTime == None:
            DeltaVar = Var - self.LastCheckTime
        else:
            DeltaVar = Var - LastCheckTime

        retVal = False

        if self.State == True:
            # Check hold time and reset state if hold time has exceeded
            # Hold indefinitely if Hold time < 0
            if self.HoldTime < 0:
                retVal = True
                return retVal
            elif (self.HoldTime + self.LastTriggerTime) <= Var:
                retVal = True # Hold to be maintained
                # Set Priority to 0 during hold
                #self.Priority = lambda x: 0
                return retVal
            else:
                self.State = False # Reset event output state
                
        if self.TriggerFlag == True:
            return False
        
        if Var >= (self.LastTriggerTime + self.HoldTime + 1000):
            self.Priority = self.PriorityOriginal

        if isinstance(self.AdditionalArguments,tuple):
            P = self.PDF(Var,DeltaVar,*self.AdditionalArguments)
        else:
            P = self.PDF(Var,DeltaVar)


        #Limit the probablity to the interval [0,1]
        P = max(min(P,1),0)

        if Estimate == True:
            # Used when only calculating the probability of occurrence from an external trigger such as information matrix factor calculation
            # Does not affect the current parameters of the event
            return P
        
        if Measure == True:
            if self.Error == False:
                return self.State
            else:
                return 0.5 # Causes InfoFactorA to compute to 0 for this event
        
        if self.CheckFixedThresholdType() == True:
            Threshold = self.Threshold(Var)
        else:
            Threshold = self.RNG.uniform()
            if self.Threshold(Var) >= 1:
                Threshold = self.Threshold(Var)
        if P >= Threshold:
            retVal = True
            self.LastTriggerTime = Var
            self.State = True
            if self.SingleTrigger == True:
                self.TriggerFlag = True
        else:
            self.State = False

        self.LastCheckTime = Var

        return retVal




if __name__ == "__main__":
    pass