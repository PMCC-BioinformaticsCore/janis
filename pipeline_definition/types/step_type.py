from abc import ABC, abstractmethod
from typing import List, Dict, Optional

from pipeline_definition.graph.node import Node, NodeType
from pipeline_definition.types.input_type import InputType, Input


class DependencySpec:
    def __init__(self, tag: Optional[str], step: Optional[str], output: Optional[str]):
        self.tag: Optional[str] = tag
        self.step: Optional[str] = step
        self.output: Optional[str] = output


class StepInput:
    def __init__(self, tag: str, input_type: str):
        """

        :param tag: tag for input, what the yml will reference (eg: input1: path/to/file)
        :param input_type:
        """
        self.tag = tag
        self.input_type = input_type


class StepOutput:
    def __init__(self, tag: str, output_type: str):
        self.tag = tag
        self.output_type = output_type


class Step(ABC):

    def __init__(self, label: str, step_meta):
        self.__id = label
        self.__meta = step_meta
        self.__inputs = {}

    # def type(self) -> InputType:
    #     return self.__type

    def id(self) -> str:
        return self.__id

    def meta(self) -> dict:
        return self.__meta

    def tool(self) -> str:
        return self.meta()["tool"]

    def identify(self):
        if True:
            print("Instance: [", self.__id, " - ", self.tool(), " - ", self.__meta, " ]")

    def input_labels(self, tag):
        return self.requires()[tag]

    def input_value(self, tag) -> str:
        return self.meta()[tag]

    @abstractmethod
    def provides(self) -> Dict[str, StepOutput]:
        # A map of output ids and output types. For example, an aligner might provide
        # {'aligned_file': bam_file}. Dependent steps can then access this as
        # align_step_id/aligned_file
        raise RuntimeError("Please provide implementation")

    @abstractmethod
    def requires(self) -> Dict[str, StepInput]:
        # Return the input types required by this step
        raise RuntimeError("Please provide implementation")

    @abstractmethod
    def translate(self, mapped_inputs) -> Dict[str, str]:
        # Return a language specific dictionary that will be translated
        # to the output text.
        raise RuntimeError("Please provide implementation")

    @abstractmethod
    def cores(self) -> int:
        # Number of CPUS
        raise RuntimeError("Please provide implementation")

    @abstractmethod
    def ram(self) -> int:
        # Amount of ram
        raise RuntimeError("Please provide implementation")

    @staticmethod
    def select_type_name_from(meta) -> str:
        return meta["tool"]

        # selection = None
        # for candidate in iter(meta.keys()):
        #     if candidate == 'tag':
        #         continue
        #     if candidate == 'input_scope':
        #         continue
        #     selection = candidate
        #     break
        # return selection

    def validate_input_output_spec(self):
        pass

        # output_specs = self.provides()
        # if output_specs:
        #   for ospec in output_specs:
        #     if not isinstance(ospec, InputType):
        #       raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")
        #
        #     if Step.STR_ID not in ospec:
        #       raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")
        #     name = ospec.get(Step.STR_ID)
        #     if not name:
        #       raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")
        #
        #     if Step.STR_TYPE not in ospec:
        #       raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")
        #     type = ospec.get(Step.STR_TYPE)
        #     if not type:
        #       raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")

    @staticmethod
    def default_step_tag_name():
        return 'untagged'

    @staticmethod
    def dependency_spec_from(requirement_value: str) -> DependencySpec:

        tag: Optional[str] = None
        step: Optional[str] = None
        output: Optional[str] = None

        parts: List[str] = requirement_value.split(".")

        head: str = parts[0]
        if head.startswith("#"):
            tag = head[1:]

        return DependencySpec(tag, step, output)
        # return {
        #     'tag': tag,
        #     'step': step,
        #     'output': output
        # }

    def provided_value_for_requirement(self, step_input: InputType):

        if self.__meta is None:
            return None

        provided = self.__meta.get(step_input.type_name())
        return provided


class StepFactory(ABC):
    @classmethod
    @abstractmethod
    def type(cls) -> str:
        pass

    @classmethod
    @abstractmethod
    def label(cls) -> str:
        pass

    @classmethod
    def description(cls) -> str:
        return cls.label()

    @classmethod
    @abstractmethod
    def schema(cls) -> dict:
        pass

    @classmethod
    @abstractmethod
    def build(cls, label: str, meta: dict) -> Step:
        pass

    @classmethod
    def support_translations(cls) -> List[str]:
        return ['cwl']

    @classmethod
    def build_from(cls, label: str, step_meta: dict, debug=False) -> Step:
        step_type = cls.type()
        if debug:
            print(step_type, "factory: Building from", step_meta)
        obj = cls.build(label, step_meta)
        obj.identify()
        return obj


class StepNode(Node):
    def __init__(self, step: Step):
        super().__init__(NodeType.TASK, step.id())
        self.step = step

    def inputs(self) -> Dict[str, StepInput]:
        return self.step.requires()

    def outputs(self) -> Dict[str, StepOutput]:
        return self.step.provides()


class TaggedDatum(ABC):
    @abstractmethod
    def tags(self):
        # A set tags that can select among similar types
        pass

    @abstractmethod
    def datum_type(self):
        # A datum_type
        pass

    def satisfies(self, datum):
        # A concrete implementation here
        pass
