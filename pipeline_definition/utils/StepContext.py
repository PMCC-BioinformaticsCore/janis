import json

class StepContext:
    def __init__(self, step):
        self.__step = step

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,sort_keys=True, indent=4)

    def print(self):
        print(self.toJSON())
