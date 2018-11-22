from examples.unix_commands import Untar
import yaml
from pipeline_definition.types.tool import ToolInput

t = Untar()

d = t.cwl()
print(yaml.dump(d))
print(t.cwl())