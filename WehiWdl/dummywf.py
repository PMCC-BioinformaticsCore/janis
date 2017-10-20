from WorkflowObjectModel import *;


class DummyWorkflow(Workflow):
    def __init__(self):

        cmd = CommandShell();
        step = Step();
        step.setTask( cmd );
        self.steps.append( step );
        pass;
