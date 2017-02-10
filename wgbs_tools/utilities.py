"""
Common utilities used to run and test wgbs_tools modules
"""

import os


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
