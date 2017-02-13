"""
Implementing the OutputTimeControler class
"""


class OutputTimeControler(object):
    """
    Compute the time or iteration of outputs
    """
    def __init__(self, identifier, time_period=None, iteration_period=None):
        self.__id = identifier
        if time_period is None and iteration_period is None:
            iteration_period = 1
        elif time_period is not None and iteration_period is not None:
            raise ValueError("Only time_period nor iteration_period can be specified for the OutputTimeControler {:s}".format(self))
        self.__time_period = time_period
        self.__iteration_period = iteration_period
        self.__next_output_time = self.__time_period
        self.__next_output_iteration = self.__iteration_period

    def __str__(self):
        return self.__class__.__name__ + " of the database {:s}".format(self.__id)

    def db_has_to_be_updated(self, time, iteration):
        """
        Return True if the iteration or time requires to write fields in the database
        and update the next output time or iteration
        """
        answer = False
        try:
            if time > self.__next_output_time:
                answer = True
                self.__next_output_time += self.__time_period
        except TypeError:
            if iteration > self.__next_output_iteration:
                answer = True
                self.__next_output_iteration += self.__iteration_period
        return answer
