from typing import Dict, List, Tuple, Set, Optional, Any

from Pipeline.types.data_types import DataType, NativeTypes
from Pipeline.utils.logger import Logger, LogLevel
from Pipeline.utils.descriptor import BaseDescriptor, GetOnlyDescriptor

import networkx as nx

from Pipeline.graph.node import Node, NodeTypes, NodeAndTag, layout_nodes2
from Pipeline.utils.errors import DuplicateLabelIdentifier, InvalidNodeIdentifier, NodeNotFound, InvalidStepsException, \
    InvalidInputsException
from Pipeline.workflow.input import Input, InputNode
from Pipeline.workflow.output import Output, OutputNode
from Pipeline.workflow.step import Step, StepNode
from Pipeline.tool.tool import Tool, ToolInput, ToolOutput


class Workflow:
    name = BaseDescriptor()

    def __init__(self, name: str):
        Logger.log(f"Creating workflow with name: '{name}'")

        self.name = name

        self._labels: Dict[str, Node] = {}

        self._inputs: List[InputNode] = []
        self._steps: List[StepNode] = []
        self._outputs: List[OutputNode] = []

        self.graph: nx.MultiDiGraph = nx.MultiDiGraph()
        self.connections: Dict[str, Set[NodeAndTag]] = {}

    def draw_graph(self):
        import matplotlib.pyplot as plt

        default_color = 'black'

        G = self.graph
        edges_attributes = [G.edges[e] for e in G.edges]
        edge_colors = [x["color"] if 'color' in x else default_color for x in edges_attributes]
        node_colors = [NodeTypes.to_col(x.node_type) for x in G.nodes]

        pos = layout_nodes2(list(G.nodes))

        for n in pos:
            G.node[n]['pos'] = pos[n]

        nx.draw(G, pos=pos, edge_color=edge_colors, node_color=node_colors, with_labels=True)
        plt.show()

    def add_input(self, inp: Input):
        Logger.log(f"Adding input '{inp.id()}' to '{self.name}'")
        node: InputNode = InputNode(inp)
        self._add_node(node)
        self._inputs.append(node)

    def add_output(self, outp: Output):
        Logger.log(f"Adding output '{outp.id()}' to '{self.name}'")
        node: OutputNode = OutputNode(outp)
        self._add_node(node)
        self._outputs.append(node)

    def add_step(self, step: Step):
        Logger.log(f"Adding step '{step.id()}' to '{self.name}'")
        node: StepNode = StepNode(step)
        self._add_node(node)
        self._steps.append(node)

    def _add_node(self, node: Node):

        if node.id() in self._labels:
            existing = self._labels[node.id()]
            message = f"Attempted to add an input with the identifier '{node.id()}' but there was already " \
                f"an {existing.node_type} node with representation: {str(existing)}"
            Logger.log(message)
            raise DuplicateLabelIdentifier(message)
        self._labels[node.id()] = node
        self.graph.add_node(node)

    def _remove_node(self, node: Node):

        if node.id() not in self._labels:
            Logger.log(f"Couldn't find node with id '{node.id()}' to remove")
            return

    @staticmethod
    def get_label_by_inference(s) -> str:
        try:
            if type(s) == str:
                return s
            else:
                return s.id()
        except Exception as e:
            Logger.log(e, LogLevel.CRITICAL)

    def add_edge(self, start, finish):

        s = self.get_label_by_inference(start)
        f = self.get_label_by_inference(finish)
        Logger.log(f"Adding edge between {s} and {f}", LogLevel.INFO)

        # set nodes
        s_parts = s.split("/")
        f_parts = f.split("/")

        # ORDER:
        #   1. Get the input type to link
        #   2. Try and find this type in the cur_step body
        #       (a) - If it's optional, skip
        #       (b) - else: throw error
        #   3. Determine where the input is coming from (input / other cur_step)
        #   4. Try to exactly connect the input to the output/tag
        #       (a) If there's only one output, log warning if it doesn't exactly match
        #       (b) If there are multiple outputs, raise Exception if it doesn't exactly match
        #   5. Check that the types match
        #   6. Add edge

        if len(s_parts) == 0 or len(s_parts) > 2:
            message = f"The start tag '{s}' contained an invalid number of parts ({len(s_parts)}), " \
                f"it should follow the format: 'label/tag' (1 or 2 parts)"
            Logger.log(message, LogLevel.CRITICAL)
            raise InvalidNodeIdentifier(message)
        if len(f_parts) == 0 or len(f_parts) > 2:
            message = f"The finish tag '{f}' contained an invalid number of parts ({len(f_parts)}), " \
                f"it should follow the format: 'label/tag' (1 or 2 parts)"
            Logger.log(message, LogLevel.CRITICAL)
            raise InvalidNodeIdentifier(message)

        s_label = s_parts[0]
        f_label = f_parts[0]

        s_node = self._labels.get(s_label)
        f_node = self._labels.get(f_label)

        if s_node is None:
            message = f"The node '{s_label}' was not found in the graph, " \
                f"please add the node to the graph before trying to add an edge"
            Logger.log(message, LogLevel.CRITICAL)
            raise NodeNotFound(message)
        if f_node is None:
            message = f"The node '{f_label}' was not found in the graph, " \
                f"please add the node to the graph before trying to add an edge"
            Logger.log(message, LogLevel.CRITICAL)
            raise NodeNotFound(message)

        # we have the two nodes

        s_node.set_depth(f_node.depth - 1)
        f_node.set_depth(s_node.depth + 1)

        s_type = self.get_tag_and_type_from_start_edge_node(s_node, s_parts, f)

        if f_node.node_type == NodeTypes.OUTPUT:
            f_node: OutputNode = f_node
            f_type = f_parts[0], s_type[1]
            f_node.output.data_type = s_type[1]
            Logger.log("Connecting to output")
        else:
            f_type = self.get_tag_and_type_from_final_edge_node(f_node, f_parts)

            if s_type is None and f_type is not None:
                s_type = self.get_tag_and_type_from_start_edge_node(s_node, s_parts, f, f_type[1])

            if s_type is None and f_type is None:
                s_type, f_type = self.guess_connection_between_nodes(s_node, f_node)

            if s_type is None or f_type is None:
                raise Exception(f"Could not identify connection for edge {s} → {f}")

        # NOW: Let's build the connection
        f_node.connection_map[f_parts[-1]] = s_type[0], s_node

        correct_type = f_type[1].can_receive_from(s_type[1])

        if not correct_type:
            Logger.log(f"Mismatch of types when joining '{s}' to {f}' "
                       f"({s_type[1].id()} -/→ {f_type[1].id()})",
                       LogLevel.CRITICAL)
            Logger.log(f"No action taken to correct type-mismatch of '{s}' "
                       f"to {f}'")

        col = 'black' if correct_type else 'r'
        self.graph.add_edge(s_node, f_node, type_match=correct_type, color=col)

    def get_tag_and_type_from_start_edge_node(self, node: Node, input_parts: List[str], referenced_by: str,
                                              guess_type: Optional[DataType] = None) -> Tuple[str, DataType]:

        if node.node_type == NodeTypes.OUTPUT:
            raise Exception(f"Can't join '{referenced_by}' to output node '{node.id()}'")

        lbl = input_parts[0]
        s = "/".join(input_parts)
        if node.node_type == NodeTypes.INPUT:
            if len(input_parts) != 1:
                message = f"The input tag '{s}' over-represents the INPUT node, this was corrected ({s} → {lbl})"
                Logger.log(message, LogLevel.WARNING)
            return lbl, next(iter(node.outputs().values())).output_type

        # We are a step node now
        node: StepNode = node
        outs: Dict[str, ToolOutput] = node.outputs()

        types = [(x, outs[x].output_type) for x in outs]

        if len(types) == 0:
            raise InvalidStepsException(f"The step '{referenced_by}' referenced the step '{node.id()}' with tool "
                                        f"'{node.step.get_tool().id()}' that has no inputs")
        elif len(types) == 1:
            if len(input_parts) != 2:
                under_over = "under" if len(input_parts) < 2 else "over"
                Logger.log(f"The node '{node.id()}' {under_over}-referenced an output of the tool "
                           f"'{node.step.get_tool().id()}', this was automatically corrected "
                           f"({s} → {lbl}/{types[0][0]})", LogLevel.WARNING)
                s = f"{lbl}/{types[0][0]}"
            return s, types[0][1]

        else:
            # if the edge tag doesn't match, we can give up
            tag = input_parts[1]
            t = node.step.get_tool().outputs_map().get(tag)
            if t:
                if len(input_parts) != 2:
                    Logger.log(f"The node '{node.id()}' did not correctly reference an output of the tool "
                               f"'{node.step.get_tool().id()}', this was automatically corrected "
                               f"({s} → {lbl}/{tag})", LogLevel.WARNING)
                    s = f"{lbl}/{tag}"
                return s, t.output_type

            possible_tags = ", ".join(f"'{x}'" for x in outs)
            raise Exception(f"Could not identify an output called '{tag}' on the node '{node.id()}' "
                            f", possible tags: {possible_tags}")

    @staticmethod
    def get_tag_and_type_from_final_edge_node(node: Node, input_parts: List[str]) -> Tuple[str, DataType]:
        if node.node_type == NodeTypes.INPUT:
            raise Exception(f"Can't connect TO input node '{node.id()}'")
        if node.node_type == NodeTypes.OUTPUT:
            raise Exception("An internal error has occurred: output nodes should be filtered from the "
                            "function 'get_tag_and_type_from_final_edge_node'")

        node: StepNode = node
        ins = node.inputs()
        types = [(x, ins[x].input_type) for x in ins]
        lbl = input_parts[0]
        s = "/".join(input_parts)

        if len(types) == 0:
            raise InvalidStepsException(f"The step '{node.id()}' has no inputs, and cannot be added to this workflow")

        if len(types) == 1:

            if len(input_parts) != 2:
                Logger.log(f"The node '{node.id()}' did not correctly reference an input of the tool "
                           f"'{node.step.get_tool().id()}', this was automatically corrected "
                           f"({s} → {lbl}/{types[0][0]})", LogLevel.WARNING)
                # s = f"{lbl}/{types[0][0]}"
            return types[0]
        elif len(input_parts) < 2:
            possible_tags = ", ".join(f"'{x}'" for x in ins)
            raise Exception(f"The tag '{s}' could not uniquely identify an input of '{node.id()}', requires the one of"
                            f"the following tags: {possible_tags}")
        else:
            tag = input_parts[1]
            t = node.step.get_tool().inputs_map().get(tag)
            if t:
                if len(input_parts) != 2:
                    Logger.log(f"The node '{node.id()}' did not correctly reference an input of the tool "
                               f"'{node.step.get_tool().id()}', this was automatically corrected "
                               f"({s} → {lbl}/{tag})", LogLevel.WARNING)
                    s = f"{lbl}/{tag}"
                return s, t.input_type

            possible_tags = ", ".join(f"'{x}'" for x in ins)
            raise Exception(f"Could not identify an input called '{tag}' on the node '{node.id()}' "
                            f", possible tags: {possible_tags}")

        raise Exception("Unhandled pathway in 'get_tag_and_type_from_node'")

    @staticmethod
    def guess_connection_between_nodes(s_node: Node, f_node: Node) -> Tuple[Tuple[str, DataType], Tuple[str, DataType]]:
        outs, ins = s_node.outputs(), f_node.inputs()

        s_types: List[Tuple[str, DataType]] = [(x, outs[x].output_type) for x in outs]
        f_types: List[Tuple[str, DataType]] = [(x, ins[x].input_type) for x in ins]

        # O(n**2) for determining types
        matching_types: List[Tuple[Tuple[str, DataType], Tuple[str, DataType]]] = []
        for s_type in s_types:
            matching_types.extend([(s_type, f_type) for f_type in f_types if f_type[1].can_receive_from(s_type)])

        if len(matching_types) == 0:
            raise InvalidInputsException(f"Can't find a common type between '{s_node.id()}' and '{f_node.id()}'")

        if len(matching_types) > 1:
            raise InvalidInputsException(f"Couldn't guess the single connection between '{s_node.id()}' "
                                         f"and '{f_node.id()}', there was more than 1 compatible connection")

        return matching_types[0]

    def attempt_connect_node(self, node: Node):
        """
        Get the inputs and outputs of the node, and try to make the connections
        :param node:
        """
        unfulfilled_connections: Dict[str, List[NodeAndTag]] = {}

        inputs: Dict[str, ToolInput] = node.inputs()
        outputs = node.outputs()

        if node.node_type == NodeTypes.INPUT:
            pass
            # try to fulfill unfulfilled connections

        if node.node_type == NodeTypes.TASK:
            pass

    def cwl(self) -> Tuple[Dict[str, Any], Dict[str, Any], List[Dict[str, Any]]]:
        # Let's try to emit CWL
        d = {
            "class": "Workflow",
            "cwlVersion": "v1.0",
            "id": self.name,
            "label": self.name,
            "requirements": [
                {"class": "InlineJavascriptRequirement"}
            ]
        }

        if self._inputs:
            d["inputs"] = {i.id(): i.cwl() for i in self._inputs}

        if self._outputs:
            d["outputs"] = {o.id(): o.cwl() for o in self._outputs}

        if self._steps:
            d["steps"] = {s.id(): s.cwl() for s in self._steps}

        tools = []
        tools_to_build: Dict[str, Tool] = {s.step.get_tool().tool(): s.step.get_tool() for s in self._steps}
        for t in tools_to_build:
            tools.append(tools_to_build[t].cwl())

        inp = {i.id(): i.input_cwl_yml() for i in self._inputs}

        return d, inp, tools

    def wdl(self):

        get_alias = lambda t: t[0] + "".join([c for c in t[1:] if c.isupper()])

        tools = [s.step.get_tool() for s in self._steps]
        tool_name_to_tool: Dict[str, Tool] = {t.tool().lower(): t for t in tools}
        tool_name_to_alias = {}
        steps_to_alias: Dict[str, str] = {s.id().lower(): get_alias(s.id()).lower() for s in self._steps}

        aliases = set()

        for tool in tool_name_to_tool:
            a = get_alias(tool).upper()
            s = a
            idx = 2
            while s in aliases:
                s = a + idx
                idx += 1
            aliases.add(s)
            tool_name_to_alias[tool] = s

        tab_char = '  '
        nline_char = '\n'

        import_str ="import \"tools/{tool_file}.wdl as {alias}"
        input_str = "{tb}{data_type} {identifier}{default_with_equals_if_required}"
        step_str =  "{tb}call {tool_file}.{tool} as {alias} {{ input: {tool_mapping} }}"
        output_str ="{tb2}{data_type} {identifier} = {alias}.{outp}"

        imports = [import_str.format(
            tb=tab_char,
            tool_file=t,
            alias=tool_name_to_alias[t.lower()].upper()
        ) for t in tool_name_to_tool]

        inputs = [input_str.format(
            tb=tab_char,
            data_type=i.input.data_type.wdl(),
            identifier=i.id(),
            default_with_equals_if_required=""
        ) for i in self._inputs]

        steps = [step_str.format(
            tb=tab_char,
            tool_file=tool_name_to_alias[s.step.get_tool().tool().lower()].upper(),
            tool=s.step.get_tool().wdl_name(),
            alias=s.id(),
            tool_mapping=', '.join(s.wdl_map()) #[2 * tab_char + w for w in s.wdl_map()])
        ) for s in self._steps]

        outputs = [output_str.format(
            tb2=2 * tab_char,
            data_type=o.output.data_type.wdl(),
            identifier=o.id(),
            alias=steps_to_alias[next(iter(o.connection_map.values()))[0].split('/')[0].lower()].lower(),
            outp=o.id()
        ) for o in self._outputs]

        # imports = '\n'.join([f"import \"tools/{t}.wdl\" as {tool_name_to_alias[t.lower()].upper()}" for t in tool_name_to_tool])
        # inputs = '\n'.join([f"{tab_char}{i.input.data_type.wdl()} {i.id()}" for i in self._inputs])
        # steps = '\n'.join([f"{tab_char}call {tool_name_to_alias[s.step.get_tool().tool().lower()].upper()}"
        #                    f".{s.step.get_tool().wdl_name()} as {s.id()} {{input: \n{(',' + nline_char).join([2 * tab_char + w for w in s.wdl_map()])}\n{tab_char}}}" for s in self._steps])
        # outputs = '\n'.join([f"{2*tab_char}{o.output.data_type.wdl()} {o.id()} = {steps_to_alias[next(iter(o.connection_map.values()))[0].split('/')[0].lower()].lower()}.{o.id()}" for o in self._outputs])

        workflow = f"""
{imports}

workflow {self.name} {{
{nline_char.join(inputs)}

{nline_char.join(steps)}

{tab_char}output {{
{nline_char.join(outputs)}
{tab_char}}}
}}"""
        tools = {t.id(): t.wdl() for t in tools}
        inp = {f"{self.name}.{i.id()}": i.input.data_type.wdl() for i in self._inputs}

        return workflow, inp, tools
