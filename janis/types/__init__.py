from abc import ABC


class Selector(ABC):
    pass


class InputSelector(Selector):

    def __init__(self, input_to_select, suffix=None, prefix=None):
        self.input_to_select = input_to_select
        self.prefix = prefix
        self.suffix = suffix


class WildcardSelector(Selector):
    def __init__(self, wildcard):
        self.wildcard = wildcard
