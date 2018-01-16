#
# Registry of objects
#

class RegistryException(Exception):
    pass

class Registry:
    def __init__(self):
        self.registry = {}

    def register(self, obj):
        if obj.type() in self.registry:
            raise RegistryException('Type %s is already registered.' % obj.type())

        self.registry[obj.type()] = obj

    def objects(self):
        return self.registry.values()

    def object(self, type_name):
        return self.registry[type_name]
