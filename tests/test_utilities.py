"""
Test suite for utilities
"""

from wgbs_tools import utilities


def test_which_returnexpected():
    """Tests which function returns expected string"""
    assert utilities.which('python'), \
        'Which function failing, Python cannot be found in path. However, ' \
        'you are running python if you are running this command, so ' \
        'something is wrong with the which command!'


def test_which_returnnone():
    """
    Tests which function to return None when prompted with something not in
    $PATH
    """
    assert not utilities.which(''), 'Which is not returning None, and it ' \
                                    'should in this test case.'


def test_show_value(sys_version):
    """Tests show_value"""
    if sys_version == 2:
        string = u'asdf'
    else:
        string = 'asdf'
    newstring = utilities.show_value(string)
    assert newstring == 'asdf', \
        'Python version string conversion for BedTools is not working ' \
        'correctly. Python major version:{}'.format(sys_version)

