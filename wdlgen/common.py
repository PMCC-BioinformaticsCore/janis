from typing import Union, List

from wdlgen.types import WdlType
from wdlgen.util import WdlBase


class Input(WdlBase):
    def __init__(self, data_type: Union[WdlType, List[WdlType]], name: str, expression: str = None):
        self.type = data_type
        self.name = name
        self.expression = expression

    def wdl(self):
        if self.type is None:
            raise Exception(f"Could not convert wdlgen.Input ('{self.name}') to string because type was null")

        f = "{type} {name}{def_w_equals}"
        if isinstance(self.type, list):

            return [
                f.format(
                    type=self.type[i].wdl(),
                    name=self.name + ("" if i == 0 else "_" + str(i)),
                    def_w_equals=(" = {val}".format(val=self.expression) if self.expression else ""))
                for i in range(len(self.type))]

        wd = self.type.wdl()
        if isinstance(wd, list):

            return [f.format(
                type=t,
                name=self.name,
                def_w_equals=(" = {val}".format(val=self.expression) if self.expression else "")
            ) for t in wd]
        return f.format(
                type=wd,
                name=self.name,
                def_w_equals=(" = {val}".format(val=self.expression) if self.expression else "")
            )


class Output(WdlBase):
    def __init__(self, data_type: WdlType, name: str, expression: str = None):
        self.type = data_type
        self.name = name
        self.expression = expression

    def wdl(self):
        f = "{type} {name}{def_w_equals}"
        if isinstance(self.type, list):
            return [
                f.format(
                    type=self.type[i].wdl(),
                    name=self.name + ("" if i == 0 else "_" + str(i)),
                    def_w_equals=(" = {val}".format(val=self.expression) if self.expression else ""))
                for i in range(len(self.type))]

        wd = self.type.wdl()
        if isinstance(wd, list):

            return [f.format(
                type=t,
                name=self.name,
                def_w_equals=(" = {val}".format(val=self.expression) if self.expression else "")
            ) for t in wd]
        return f.format(
                type=wd,
                name=self.name,
                def_w_equals=(" = {val}".format(val=self.expression) if self.expression else "")
            )
