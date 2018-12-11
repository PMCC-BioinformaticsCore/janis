from typing import Dict, Any

from Pipeline.types.common_data_types import String
from Pipeline.types.data_types import NativeType, NativeTypes


class Filename(String):

    def __init__(self, extension: str=None):
        """
        :param extension: with no '.' (dot)
        """
        super().__init__(optional=True)
        self.extension = extension

    @staticmethod
    def name() -> str:
        return "Filename"

    @staticmethod
    def primitive() -> NativeType:
        return NativeTypes.kStr

    @staticmethod
    def doc() -> str:
        return """
This class is a placeholder for generated filenames, by default it is optional and CAN be overrided, however the program has been structured in a way such that these names will be generated based on the step label. These should only be used when the tool _requires_ a filename to output and you aren't concerned what the filename should be.
The Filename DataType should NOT be used as an output.
""".strip()

    @classmethod
    def schema(cls) -> Dict:
        pass

    def extension_if_required_with_dot(self) -> str:
        return ("." + self.extension) if self.extension is not None else ""

    def cwl(self) -> Dict[str, Any]:
        return {
            **super().cwl(),
            "default": self.generated_filename()
        }

    def generated_filename(self) -> str:
        import uuid
        return "filename-generated-" + str(uuid.uuid1()) + self.extension_if_required_with_dot()

    def default(self) -> str:
        return self.generated_filename()