#########################################################################
#                                                                       #
#        Python implementation of a double crystal ball function,       #
#           its integral, and the derivative of it's logarithm          #
#                       Kaustuv Datta; 08/2023                          #
#                                                                       #
#########################################################################
#                                                                       #
#                     Original C++ implementation in                    # 
# dasanalysissystem.docs.cern.ch/namespaceDAS_1_1DoubleCrystalBall.html #
#                                                                       #
#########################################################################

import math
import sys
import ROOT
from ROOT import *
from array import array

# Constants 
dinf = float('inf')
deps = sys.float_info.epsilon

class Base:
    # Constructor
    def __init__(self, f=None):
        
        # basic parameters
        self.N = 1. # global normalisation (best is to fix it to one)
        self.mu = 0. # mean
        self.sigma = 1. # width
        self.kL = float('-inf') # LHS turn-on poin, -infinity in Python
        self.nL = 2. # LHS power
        self.kR = float('inf')  # RHS turn-on point, +infinity in Python
        self.nR = 2. # RHS power
        
        # normalised parameters
        self.z = 0. # normalised resolution
        self.aR = float('inf') # normalised RHS turn-on point
        self.aL = float('-inf') # normalised LHS turn-on point
        

        # values of core at transition points
        self.expaR2 = 0. # value of core at normalised RHS turn-on point
        self.expaL2 = 0. # value of core at normalised LHS turn-on point
        
        # normalisation
        self.CR = 0. # area in RHS tail
        self.CL = 0. # area in LHS tail
        self.DD = math.sqrt(2. * math.pi) # area in core
        
        
        # derived params
        self.result = 0. # stores last evaluation
        #######################

        self.AR = float('inf')  # RHS tail normalisation
        self.AL = float('inf') # LHS tail normalisation
        self.BR = float('inf')  # RHS offset
        self.BL = float('inf') # LHS offset

        # these parameters are usually used in the analytical formula
        # but in practice the n^n factor can get too large 
        # therefore we actually do not use them in practice
        # but we still calculate them 
        
        #######################

        self.K = 0. # total normalization


        if f is not None:
            if len(f) != 7:
                raise AssertionError("f should have 7 parameters")

            # If a function is provided, then the parameters are taken from f
            p = f.GetParameters()
            self.SetParameters(p)

        
    def dump(self):
        # Code for dumping internal parameters
        # Just printing all internal parameters AND the last result
        
        print(self.N, self.mu, self.sigma, self.kL, self.nL, self.kR, self.nR,
              self.z, self.aR, self.aL, 
              self.expaR2, self.expaL2, 
              self.CR, self.CL, self.DD, self.result)

    def Eval(self, x):
        raise NotImplementedError("Eval method should be implemented by the derived class")
        #pass  # Implemented by derived classes

    def SetParameters(self, p):
        # Code for setting parameters

        self.N = p[0]
        self.mu = p[1]
        self.sigma = max(1e-8, abs(p[2]))
        assert self.sigma > 0
        self.kL = min(self.mu - deps, p[3])
        self.nL = max(1 + deps, p[4])  # Assuming `deps` in the C++ code is 1e-8
        self.kR = max(self.mu + 1e-8, p[5])
        self.nR = max(1 + deps, p[6])

        if self.aL < deps:
            self.aL = deps
        elif self.aR < deps:
            self.aR = deps
            
        self.aR = abs(self.Z(self.kR))
        self.aL = abs(self.Z(self.kL))

        self.expaR2 = math.exp(-pow(self.aR, 2) / 2.)
        self.expaL2 = math.exp(-pow(self.aL, 2) / 2.)

        self.AR = pow(self.nR / self.aR, self.nR) * self.expaR2
        self.AL = pow(self.nL / self.aL, self.nL) * self.expaL2
        self.BR = self.nR / self.aR + self.aR
        self.BL = self.nL / self.aL + self.aL

        self.CR = self.nR * self.expaR2 / (self.aR * (self.nR - 1))
        self.CL = self.nL * self.expaL2 / (self.aL * (self.nL - 1))
        self.DD = math.sqrt(math.pi / 2.) * (math.erf(self.aR / math.sqrt(2)) + math.erf(self.aL / math.sqrt(2)))
        self.K = self.CR + self.CL + self.DD

    def Z(self, x):
        return (x - self.mu) / self.sigma

    def operator(self, x, p):
        self.SetParameters(p)
        self.z = self.Z(x[0])  
        self.Eval(x)
        return self.result
    def NDim(self):
        return 1

    def wrapped_eval(self, x, p):
        self.SetParameters(p)
        return self.operator(x, p)
    
    def Clone(self):
        clone = self.__class__()  # Generalization to create an instance of the current class
        ROOT.SetOwnership(clone, False)
        return clone

    def DoEval(self, x):
        self.Eval(x)
        return self.result
    
class Distribution(Base, ROOT.Math.IMultiGenFunction):
    
    def __init__(self, f=None):
        Base.__init__(self, f)
        #ROOT.Math.IMultiGenFunction.__init__(self)

    def Eval(self, x):
        x_val = x[0]  
        if x_val < self.kL:
            R = self.aL / self.nL
            #print(f"R: {R}, self.aL: {self.aL}, self.z: {self.z}, self.nL: {self.nL}")
            #print(f"Inside power term: {1 - R * (self.aL + self.z)}")
            self.result = (self.N / self.K) * self.expaL2 * math.pow(1 - R * (self.aL + self.z), -self.nL) / self.sigma
        elif x_val < self.kR:
            self.result = (self.N / self.K) * math.exp(-math.pow(self.z, 2) / 2.) / self.sigma
        else:
            R = self.aR / self.nR
            self.result = (self.N / self.K) * self.expaR2 * math.pow(1 - R * (self.aR - self.z), -self.nR) / self.sigma
    
    
class DLog(Base, ROOT.Math.IMultiGenFunction):
    
    def __init__(self, f=None):
        Base.__init__(self, f)

    def Eval(self, x):
        x_val = x[0]  
        if x_val < self.kL:
            B = self.nL / self.aL - self.aL
            self.result = (self.nL / self.sigma) / (B - self.z)
        elif x_val < self.kR:
            self.result = -self.z / self.sigma
        else:
            B = self.nR / self.aR - self.aR
            self.result = -(self.nR / self.sigma) / (B + self.z)
    
    
class Integral(Base, ROOT.Math.IMultiGenFunction):
    
    def __init__(self, f=None):
        Base.__init__(self, f)

    def Eval(self, x):
        x_val = x[0]  
        if x_val < self.kL:
            R = self.aL / self.nL
            self.result = (self.N / self.K) * self.expaL2 * math.pow(1 - R * (self.aL + self.z), 1 - self.nL) / (R * (self.nL - 1))
        elif x_val < self.kR:
            self.result = (self.N / self.K) * (self.CL + math.sqrt(math.pi / 2.) * (math.erf(self.z / math.sqrt(2)) + math.erf(self.aL / math.sqrt(2))))
        else:
            R = self.aR / self.nR
            self.result = (self.N / self.K) * (self.CL + self.DD + self.expaR2 * (1. - math.pow(1 - R * (self.aR - self.z), 1 - self.nR)) / (R * (self.nR - 1)))

    