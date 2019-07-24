from janis import Workflow, Input, Step, Output, String
from janis.unix.tools.echo import Echo


class EdgesPreferred(Workflow):
    """
    METHOD 1 (Preferred)

    This workflow represents the recommended way to connect
    inputs, steps and outputs together, through directly
    referencing the node, or using the dot notation.

    This method removes the need to call the Workflow.add_items
    method, and is also resilient to ID changes of the nodes.
    The program will also give you more immediate feedback if
    the tag on the step is not correct (eg: if you mispelled 'inp').

    """

    def __init__(self):
        super().__init__("edgeWorkflowExample")

        inp = Input("inputId", data_type=String())
        stp = Step("stepId", tool=Echo())
        out = Output("outputId")

        self.add_edge(inp, stp.inp)
        self.add_edges([(stp.out, out)])


# Method 2 - Prone to issues when IDs change


class EdgesSecondary(Workflow):
    """
    METHOD 2 (Not recommended)

    This workflow shows a secondary method to create edges
    where directly accessing the physical input / stp / out
    may not be possible or practical.

    This method is not recommended as its not resilient to
    ID changes.
    """

    def __init__(self):
        super().__init__("edgeWorkflowExample")

        inp = Input("inputId", data_type=String())
        stp = Step("stepId", tool=Echo())
        out = Output("outputId")

        self.add_items(inp, stp, out)
        self.add_edge("inputId", "stepId/inp")
        self.add_edges([("stepId/out", "outputId")])


if __name__ == "__main__":
    w1, _, _ = EdgesPreferred().translate("wdl", to_console=False)
    w2, _, _ = EdgesSecondary().translate("wdl", to_console=False)

    print("Same output: " + str(w1 == w2))
