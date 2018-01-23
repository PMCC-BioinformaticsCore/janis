import argparse
import atexit

import pipeline_definition.bio_informatics

from pipeline_definition.pdx import PDX


def main( opts ):
    def atExit():
        print("BYE!!")
    atexit.register( atExit )

    #Get specified file to translate
    pdfile = opts.pdfile

    pdfile = "pd_1.yml"

    pdx = PDX()
    #pdx.translate( pdfile, outfile="/Users/mohammadbhuyan/Temp/out.pdx", overwriteOutfile=True  )
    pdx.translate(pdfile, outfile="/Users/mohammadbhuyan/Temp/out.pdx")


if __name__ == "__main__":
    argprsr = argparse.ArgumentParser()
    argprsr.add_argument('pdfile', help='Pipeline Definition file.')
    opts = argprsr.parse_args()
    main( opts );

