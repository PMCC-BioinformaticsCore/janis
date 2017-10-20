import sys;
import os;

class InputOutputSchema():
    def __init__(self):
        pass;

class Task():
    def __init__(self):
        pass;
    def emitCWL(self, file, logger):
        file.write("X");
        pass;

class CommandShell(Task):
    def __init__(self):
        pass;

class Step():
    _task = None;
    def __init__(self):
        pass;

    def setTask(self, task):
        self._task = task;

    def emitCWL(self, file, logger):
        self._task.emitCWL( file, logger );
        pass;

class Documentation():
    pass;


class Workflow():
    steps = [];

    def __init__(self):
        pass;

    def emitCWL(self, file, logger ):
        with open( file, 'w+' ) as fh:
            logger.info("CWL FILE:" + os.path.abspath(file));

            for step in self.steps:
                step.emitCWL( fh, logger );
            pass;
        pass;

