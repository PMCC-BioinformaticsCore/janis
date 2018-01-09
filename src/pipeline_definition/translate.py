
import sys
import yaml


def load_input():
    if len(sys.argv) < 1:
        sys.exit('An input file name must be specified.')

    fn = sys.argv[1]
    with open(fn, 'r') as f:
        data = f.read()

    return data


def main():
    input_def = load_input()
    yml = yaml.load(input_def)
    print(yml)


if __name__ == '__main__':
    main()

