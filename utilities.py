"""
Common utilities used to run and test wgbs_tools modules
"""

import os
from collections import defaultdict
import sys

def which(program):
    """
    A python version of the linux `which` command.

    Useful bits of function copied from:
    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    :param program: Name of program that should be accessible in path
    :return: Path to program or None if not in $PATH
    """
    def is_exe(fpath):
        """Small function that determines if something is executable"""
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath = os.path.split(program)[0]
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def find_occurences(string, char):
    """
    Copied from: http://stackoverflow.com/questions/13009675/
        find-all-the-occurrences-of-a-character-in-a-string

    :param string: string to be searched
    :param char: character to search
    :return: list of positions of ocurrences of character in string
    """
    return [it for it, letter in enumerate(string) if letter == char]


def nested_dict(n, type):
    """
    Creates a nested dictionary.
    Copied from: http://stackoverflow.com/questions/29348345/
        declaring-a-multi-dimensional-dictionary-in-python

    :param n: number of nested levels
    :param type: type of data, ex: str
    :return: returns the dict
    """
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

def show_value(string):
    """
    Convert unicode to str under Python 2; all other values pass through
    unchanged. This is necessary for pybedtools. Code found at:
    https://daler.github.io/pybedtools/intervals.html

    :param string: str that needs to be converted
    :return: converted str
    """
    if sys.version_info.major == 2:
        if isinstance(string, unicode):
            return str(string)
    return string


