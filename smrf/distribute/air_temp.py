__author__ = "Scott Havens"
__maintainer__ = "Scott Havens"
__email__ = "scott.havens@ars.usda.gov"
__date__ = "2015-12-30"

# import numpy as np
import logging
from smrf.distribute import image_data

#import matplotlib.pyplot as plt

class ta(image_data.image_data):
    """
    The ``ta()`` class allows for variable specific distributions that 
    go beyond the base class.
    
    Air temperature 
    
    Attributes:
        config: configuration from [air_temp] section
        air_temp: numpy matrix of the air temperature
        stations: stations to be used in alphabetical order
    
    """
    
    variable = 'air_temp'
    
    # these are variables that can be output
    output_variables = {'air_temp': {
                                     'units': 'degree Celcius',
                                     'long_name': 'air_temperature'
                                     }
                        } 
    
    def __init__(self, taConfig):
        
        # extend the base class
        image_data.image_data.__init__(self, self.variable)
        self._logger = logging.getLogger(__name__)
        
        # check and assign the configuration
        self.getConfig(taConfig)
        
        self._logger.debug('Created distribute.air_temp')
        
        
    def initialize(self, topo, metadata):
        """
        Initialize the distribution, calls image_data.image_data._initialize()
        
        Args:
            topo: smrf.data.loadTopo.topo instance contain topo data/info
            metadata: metadata dataframe containing the station metadata
                        
        """
        
        self._logger.debug('Initializing distribute.air_temp')
        self._initialize(topo, metadata)
        
                     
        
    def distribute(self, data):
        """
        Distribute air temperature
        
        Args:
            data: air_temp data frame
            
        """
    
        self._logger.debug('%s -- Distributing air_temp' % data.name)
        
        self._distribute(data)
        
        
    def distribute_thread(self, queue, data):
        """
        Distribute the data using threading and queue
         
        Args:
            queue: queue dict for all variables
            data: pandas dataframe for all data required
         
        Output:
            Changes the queue air_temp for the given date
        """
         
        for t in data.index:
             
            self.distribute(data.ix[t])
         
            queue[self.variable].put( [t, self.air_temp] )
    

    
    
    