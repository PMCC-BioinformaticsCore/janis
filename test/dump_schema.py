#
# Test loading the registry and dumping the schema as YAML
#
import sys
from pipeline_definition import translate

fn = sys.argv[1]

translate(fn)