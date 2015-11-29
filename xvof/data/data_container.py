'''
Created on 29 nov. 2015

@author: guillaume2
'''
import lxml.etree as et

class Singleton(type):
    """
    A metaclass for creating singleton object
    """
    instance = None
    def __call__(cls, *args, **kw):
        if not cls.instance:
            cls.instance = super(Singleton, cls).__call__(*args, **kw)
        return cls.instance

class DataContainer(object):
    """
    Contains the data read from the datafile
    """
    __metaclass__ = Singleton
    #
    def __init__(self, datafile_name="XDATA.xml"):
        '''
        Constructor
        '''
        self.__datafile = datafile_name
        self.__datadoc = et.parse(self.__datafile)
    
    def getFinalTime(self):
        return float(self.__datadoc.find('time-management/final-time').text)
    
    def getInitialTimeStep(self):
        return float(self.__datadoc.find('time-management/initial-time-step').text)
    
    def getLength(self):
        return float(self.__datadoc.find('geometry/length').text)

    def getSection(self):
        return float(self.__datadoc.find('geometry/section').text)
    
    def getEquationOfState(self):
        return self.__datadoc.find('matter/equation-of-state').text
    
    def getInitialPressure(self):
        return float(self.__datadoc.find('matter/initialization/initial-pressure').text)

    def getInitialTemperature(self):
        return float(self.__datadoc.find('matter/initialization/initial-temperature').text)

    def getInitialDensity(self):
        return float(self.__datadoc.find('matter/initialization/initial-density').text)

    def getInitialInternalEnergy(self):
        return float(self.__datadoc.find('matter/initialization/initial-internal-energy').text)
    
    def getNumberOfElements(self):
        return float(self.__datadoc.find('numeric-parameters/number-of-elements').text)

    def getLinearPseudoParameter(self):
        return float(self.__datadoc.find('numeric-parameters/linear-pseudo').text)

    def getQuadraticPseudoParameter(self):
        return float(self.__datadoc.find('numeric-parameters/quadratic-pseudo').text)
    
    def getCFL(self):
        return float(self.__datadoc.find('numeric-parameters/cfl').text)
    
    def hasExternalSolver(self):
        if self.__datadoc.find('numeric-parameters/external-solver-library') is not None:
            return True
    
    def getExternalSolverPath(self):
        return self.__datadoc.find('numeric-parameters/external-solver-library').text
    
    def getNumerOfImages(self):
        return int(self.__datadoc.find('output/number-of-images').text)
    
    def hasImagesDump(self):
        if (self.__datadoc.find('output/dump-images').text).lower() == 'true':
            return True
        else:
            return False
    
    def hasImagesShow(self):
        if (self.__datadoc.find('output/show-images').text).lower() == 'true':
            return True
        else:
            return False