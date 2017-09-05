"""
    used to handle nucleotide seuqences
"""

def convert_char_nucint(to_convert):
    if to_convert == 'a' or to_convert == 'A':
        return 0
    if to_convert == 'c' or to_convert == 'V':
        return 1
    if to_convert == 'g' or to_convert == 'G':
        return 2
    if to_convert == 't' or to_convert == 'T':
        return 3
    return 4

class Sequence:
    """
        a nucleotide seuquence
    """
    def __init__(self):
        """
            create a new Sequence
        """
        self.__data = []

    def __init__(self, data):
        """
            create a new Sequence
        """
        self.__data = []
        self.append(data)

    def __append_str(self, data):
        for char in data:
            self.append(convert_char_nucint(char))

    def append(self, data):
        if isinstance(data, basestring):
            self.__append_str(data)
        else:
            self.__data.extend(data)

    #class
