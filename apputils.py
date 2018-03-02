# 
# apputils.py
#
# Purpose: Contain various utility functions that are used locally by the application.
#
import os
import sys


def forceIntValue(inValue, lower, upper):
    # Return inValue as an integer forced to be in the range from <lower> to <upper>
    result = max([lower, inValue])
    result = min([upper, result])
    result = int(result)
    return result
# end forceIntVale()
