from typing import Dict

from janis_core import DynamicWorkflow, Array, String
from janis_unix.tools import Echo


class MyFirstDynamicWorkflow(DynamicWorkflow):
    def friendly_name(self):
        return "My first dynamic workflow"

    def id(self) -> str:
        return "MyFirstDynamicWorkflow"

    def constructor(self, inputs: Dict[str, any], hints: Dict[str, str]):

        my_value_to_print = inputs.get("inp")

        if isinstance(my_value_to_print, list):
            # print each element if 'inp' is a list
            self.input("inp_list", Array(String))
            self.step("print", Echo(inp=self.inp_list), scatter="inp")
        else:
            # print only once
            self.input("inp", String)
            self.step("print", Echo(inp=self.inp))

        # Janis will derive the output from the source (Array(String) | String)
        self.output("out", source=self.print)

    def modify_inputs(self, inputs, hints: Dict[str, str]):
        my_value_to_print = inputs.get("inp")
        if isinstance(my_value_to_print, list):
            return {**inputs, "inp_list": my_value_to_print}
        else:
            return inputs
