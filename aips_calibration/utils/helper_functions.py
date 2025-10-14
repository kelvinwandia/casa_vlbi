
import time
import logging
import os, glob, subprocess
from casatasks import *
from casaplotms import *
import bdsf
import casalogger
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def time_execution(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        
        if execution_time < 60:  # If execution time is less than a minute
            time_unit = "seconds"
            formatted_time = execution_time
        elif execution_time < 3600:  # If execution time is less than an hour
            time_unit = "minutes"
            formatted_time = execution_time / 60
        else:  # If execution time is an hour or more
            time_unit = "hours"
            formatted_time = execution_time / 3600
            
        logging.info(f"======>>>EXECUTION TIME for {func.__name__}: {formatted_time:.2f} {time_unit}")
        return result
    return wrapper

