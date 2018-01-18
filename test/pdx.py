import argparse

import os
import logging
import json
import uuid
import sys
import atexit
from signal import *
import time
import yaml
from cerberus import Validator

class AppException(Exception):
    pass

class PDX:

    def __dumpYaml(self, doc ):
        # Diagnostic - what have we got?
        print("YAML DOC [\n")
        print(yaml.dump(doc))
        print("\n] END YAML DOC \n")

    def validateSchema(self, yamlDoc ):
        self.__dumpYaml( yamlDoc )



        pass;

    def translate(self, pdfile ):
        pdfilePath = os.path.abspath(pdfile)
        print("Using PD file: " + pdfilePath)

        if not os.path.isfile(pdfilePath):
            raise ValueError("PD file does not exists.")

        # Open file stream
        file = open(pdfile, 'r')

        # Validate YAML syntax
        doc = yaml.load(file)

        # Do schema validation
        self.validateSchema( doc )


def main( opts ):
    def atExit():
        print("BYE!!")
    atexit.register( atExit )

    #Get specified file to translate
    pdfile = opts.pdfile

    pdx = PDX()
    pdx.translate( pdfile )


if __name__ == "__main__":
    argprsr = argparse.ArgumentParser()
    argprsr.add_argument('pdfile', help='Pipeline Definition file.')
    opts = argprsr.parse_args()
    main( opts );


#Python exceptions
#https://docs.python.org/2/library/exceptions.html#exception-hierarchy
