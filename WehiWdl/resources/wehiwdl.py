import sys

def identify():

    print("WDL");

class WdlObj:
    def __init__(self):
        pass;

    def identify(self):
        print(self.__class__);

class Meta(WdlObj):
    def __init__(self):
        pass;

class Resource(WdlObj):
    def __init__(self):
        pass;

class InputSpec(Resource):
    def __init__(self):
        pass;

class TaskInputSpec(InputSpec):
    def __init__(self, id, meta ):
        pass;

class TaskInputSet(WdlObj):
    def __init__(self):
        pass;
    def addInputSpec(self, inputSpec ):
        pass;

class OutputSpec(Resource):
    def __init__(self):
        pass;

class TaskOutputSpec(OutputSpec):
    def __init__(self, id, meta ):
        pass;

class TaskOutputSet(WdlObj):
    def __init__(self):
        pass;
    def addOutputSpec(self, outputSpec ):
        pass;

class Command():
    _template = None;
    def __init__(self, template):
        pass;

class TaskDocumentation(Meta):
    def __init__(self, desc, author):
        pass;

class Task(Resource):
    _documentation = None;
    def __init__(self):
        pass;

    def setDocumentation(self, doc):
        self._documentation = doc;

    def setInputSet(self, set):
        pass;

    def setOutputSet(self, set):
        pass;

    def setCommand(self, cmd ):
        pass;

class WorkflowDocumentation(Meta):
    def __init__(self, desc, author):
        pass;

class WorkflowInputSpec(InputSpec):
    def __init__(self, id, meta ):
        pass;

class WorkflowInputSet(WdlObj):
    def __init__(self):
        pass;
    def addInputSpec(self, inputSpec ):
        pass;

class WorkflowOutputSpec(OutputSpec):
    def __init__(self, id, meta ):
        pass;

class WorkflowOutputSet(WdlObj):
    def __init__(self):
        pass;
    def addOutputSpec(self, outputSpec ):
        pass;

class Workflow(WdlObj):
    _documentation = None;
    def __init__(self, label ):
        pass;

    def setDocumentation(self, doc):
        self._documentation = doc;

    def addTask(self, step, task):
        pass;

    def setInputSet(self, set):
        pass;

    def setOutputSet(self, set):
        pass;

    def addWorkflow(self, step, workflow):
        pass;

    def generateDefinition(self, format, meta ):
        pass;