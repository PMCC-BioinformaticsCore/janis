import argparse
import atexit

# Import library of already defined inouts, steps and outputs
import pipeline_definition.bio_informatics

from pipeline_definition.pdx import PDX

#Extend the translation system with user defined
from test.user_defined import UserDefinedStepFactory
from pipeline_definition.types.type_registry import register_step_factory
register_step_factory(UserDefinedStepFactory())


def main( opts ):
    def atExit():
        print("BYE!!")
    atexit.register( atExit )

    #Get specified file to translate
    pdfile = opts.pdfile

    pdfile = "pd_1.yml"

    pdx = PDX()
    pdx.translate( pdfile, outfile="/Users/mohammadbhuyan/Temp/out.pdx", overwriteOutfile=True )
    #pdx.translate(pdfile, outfile="/Users/mohammadbhuyan/Temp/out.pdx")


if __name__ == "__main__":
    argprsr = argparse.ArgumentParser()
    argprsr.add_argument('pdfile', help='Pipeline Definition file.')
    opts = argprsr.parse_args()
    main( opts )

