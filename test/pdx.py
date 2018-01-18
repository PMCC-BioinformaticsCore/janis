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
import json

from pipeline_definition.types.schema import schema

class PdxException(Exception):
    pass

class PDX:

    def __dumpYaml(self, doc ):
        # Diagnostic - what have we got?
        print("YAML DOC [\n")
        print(yaml.dump(doc))
        print("\n] END YAML DOC \n")

    def __dumpSchema(self, schema ):
        print("PDX SCHEMA [\n")
        #print(schema)
        print(json.dumps(schema, indent=4))
        print("\n] END PDX SCHEMA\n")


    def validateSchema(self, yamlDoc ):
        self.__dumpYaml( yamlDoc )

        sch = schema();
        self.__dumpSchema( sch )

        v = Validator(sch);

        v.validate(yamlDoc);


    def translate(self, pdfile, outfile=None ):
        pdfilePath = os.path.abspath(pdfile)
        print("Using PD file: " + pdfilePath)

        if not os.path.isfile(pdfilePath):
            raise ValueError("PD file does not exists.")

        if outfile is None:
            outfile = pdfile + ".pdx"
        outfilePath = os.path.abspath(outfile)

        print("Using Output file: " + outfilePath)
        if os.path.isfile(outfilePath) :
            raise ValueError("Outfile already exists.")

        if os.path.isdir(outfilePath) :
            raise ValueError("Directoiry already exists with the name of outfile.")


        # Open file stream
        file = open(pdfile, 'r')

        # Validate YAML syntax
        doc = yaml.load(file)

        # Do schema validation
        self.validateSchema( doc )


