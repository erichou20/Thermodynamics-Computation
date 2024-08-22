"""
Author: Eric Hou
Version: 1.0
"""

def lee_Kesler(M0, M1, w):
    """
    This function uses the Lee-Kesler Correlation to calculate various properties
    given the tabulated values
    :param M0: Value of the property 0
    :param M1: Value of the property 1
    :param w: Acentric factor
    :return: Estimated Value of the property
    """
    return M0 + w * M1