from janis import Array, ToolInput, ToolOutput, WildcardSelector, File
from janis.unix.data_types.tarfile import TarFile
from janis.unix.tools.unixtool import UnixTool


class Untar(UnixTool):
    @staticmethod
    def tool():
        return "untar"

    def friendly_name(self):
        return "Untar archive"

    @staticmethod
    def base_command():
        return ["tar", "xf"]

    def inputs(self):
        return [ToolInput("tarFile", TarFile(), position=0)]

    def outputs(self):
        return [ToolOutput("out", Array(File()), glob=WildcardSelector("*.java"))]
