#
# Registry of all input types
#

__registry = {}


class RegistryException(Exception):
    pass


def register_factory(factory):
    if factory.name in __registry:
        raise RegistryException('Type {name} is already registered.' % factory.name)

    __registry[factory.name] = factory


def get_factory(name):
    return __registry[name]


