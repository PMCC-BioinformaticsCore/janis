from wdlgen.types import WdlType
from wdlgen.util import WdlBase


class Input(WdlBase):
    def __init__(self, data_type: WdlType, name: str, expression: str=None):
        self.type = data_type
        self.name = name
        self.expression = expression

    def wdl(self):
        if self.type is None:
            raise Exception(f"Could not convert wdlgen.Input ('{self.name}') to string because type was null")
        return "{type} {name}{def_w_equals}".format(
            type=self.type.wdl(),
            name=self.name,
            def_w_equals=(" = {val}".format(val=self.expression) if self.expression else "")
        )


class Output(WdlBase):
    def __init__(self, data_type: WdlType, name: str, expression: str=None):
        self.type = data_type
        self.name = name
        self.expression = expression

    def wdl(self):
        return "{type} {name}{def_w_equals}".format(
            type=self.type.wdl(),
            name=self.name,
            def_w_equals=(" = {val}".format(val=self.expression) if self.expression else "")
        )
