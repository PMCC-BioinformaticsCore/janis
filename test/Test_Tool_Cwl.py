from examples.unix_commands import Untar
import yaml

t = Untar()

d = t.cwl()
print(yaml.dump(d))
print(t.cwl())