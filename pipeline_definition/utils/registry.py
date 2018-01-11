#
#
#


class RegistryException(Exception):
    pass


class Registry:
    def __init__(self):
        self.registry = {}

    def register_factory(self, factory):
        if factory.type() in self.registry:
            raise RegistryException('Type %s is already registered.' % factory.type())

        self.registry[factory.type()] = factory

    def factories(self):
        return self.registry.values()

    def factory(self, type_name):
        return self.factory[type_name]
