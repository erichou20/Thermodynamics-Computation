"""
Author: Eric Hou
Version: 1.0
"""


def standard_Prop(vals, vs):
    """
    This function takes in array or list of property values and there respective
    stoichiometric coefficients representing the species of reaction to calculate
    the overall property of the reaction
    :param vals: array of properties
    :param vs: array of stoichiometric coefficients
    :return: overall property of the reaction
    """
    if len(vals) != len(vs):
        raise IndexError("Values and Coefficients must be the same length")
    sum = 0
    for i in range(len(vals)):
        sum += vals[i] * vs[i]
    return sum
