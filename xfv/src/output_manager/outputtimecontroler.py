"""
Implementing the OutputTimeControler class
"""


class OutputTimeControler:
    """
    Compute the time or iteration of outputs
    """
    def __init__(self, identifier, time_period=None, iteration_period=None):
        self.__id = identifier
        if time_period is None and iteration_period is None:
            iteration_period = 1
        elif time_period is not None and iteration_period is not None:
            message = "Only time_period nor iteration_period can be specified "
            message += "for the OutputTimeControler {:s}".format(self)
            raise ValueError(message)

        self.__time_period = time_period
        self.__iteration_period = iteration_period
        self.__next_output_time = None
        self.__next_output_iteration = None
        if self.__time_period is not None:
            self.__next_output_time = self.__time_period
        if self.__iteration_period is not None:
            self.__next_output_iteration = self.__iteration_period

    def __str__(self):
        return (self.__class__.__name__
                + " of the database {:s} with iteration period : {} or time period {}"
                .format(self.__id, self.__iteration_period, self.__time_period))

    def db_has_to_be_updated(self, time, iteration):
        """
        Return True if the iteration or time requires to write fields in the database
        and update the next output time or iteration
        """
        answer = False
        if self.__next_output_time is not None and time >= self.__next_output_time:
            answer = True
            self.__next_output_time += self.__time_period
        if self.__next_output_iteration is not None and iteration >= self.__next_output_iteration:
            answer = True
            self.__next_output_iteration += self.__iteration_period
        return answer if time != 0. else True
