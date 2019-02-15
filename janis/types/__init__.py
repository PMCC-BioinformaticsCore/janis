from abc import ABC


class Selector(ABC):
    pass


class InputSelector(Selector):

    def __init__(self, input_to_select, extra_text=None):
        self.input_to_select = input_to_select
        self.extra_text = extra_text


class WildcardSelector(Selector):
    def __init__(self, wildcard):
        self.wildcard = wildcard
