import argparse

import atexit

from pipeline_definition.types.input_type import InputFactory
from pdx import PDX

from pipeline_definition.types.type_registry import register_input_factory

class PairedReadFactory(InputFactory):
    @classmethod
    def describe(cls):
        return {
            'forward-pattern':  {'type': 'string'},
            'backward-pattern': {'type': 'string'}
        }

    @classmethod
    def build(cls, yml):
        return None

    @classmethod
    def type(cls):
        return 'SequenceReadArchivePaired'

    @classmethod
    def label(cls):
        return 'paired read files'

    @classmethod
    def description(cls):
        return cls.label()


class BAMFactory(InputFactory):
    @classmethod
    def describe(cls):
        return {
            'path':  {'type': 'string'}
        }

    @classmethod
    def build(cls, yml):
        return None

    @classmethod
    def type(cls):
        return 'BAM'

    @classmethod
    def label(cls):
        return 'BAM file'

    @classmethod
    def description(cls):
        return cls.label()



def main( opts ):
    def atExit():
        print("BYE!!")
    atexit.register( atExit )

    #Get specified file to translate
    pdfile = opts.pdfile

    # Extend translator by registering factories
    register_input_factory(PairedReadFactory())
    register_input_factory(BAMFactory())

    pdx = PDX()
    pdx.translate( pdfile )


if __name__ == "__main__":
    argprsr = argparse.ArgumentParser()
    argprsr.add_argument('pdfile', help='Pipeline Definition file.')
    opts = argprsr.parse_args()
    main( opts );

