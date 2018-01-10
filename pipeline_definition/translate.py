
import yaml


class InvalidInput(Exception):
    pass


def translate(fn):

    with open(fn, 'r') as f:
        data = f.read()

    yml = yaml.load(data)
    build_input(yml)


def build_input(yml):
    try:
        inputs = yml['inputs']
    except KeyError:
        raise InvalidInput('The workflow must specify inputs')

    if len(inputs) == 0:
        raise InvalidInput('The workflow must specify inputs')

    print(yml)
