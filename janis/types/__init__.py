from abc import ABC


class Selector(ABC):
    pass


class InputSelector(Selector):

    def __init__(self, input_to_select, suffix=None, prefix=None, use_basename=None):
        self.input_to_select = input_to_select
        self.prefix = prefix
        self.suffix = suffix
        self.use_basename = use_basename


class WildcardSelector(Selector):
    def __init__(self, wildcard):
        self.wildcard = wildcard


class MemorySelector(Selector):
    def __init__(self, suffix=None, prefix=None):
        self.suffix = suffix
        self.prefix = prefix


class CpuSelector(Selector):
    pass
