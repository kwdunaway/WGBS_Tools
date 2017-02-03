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
