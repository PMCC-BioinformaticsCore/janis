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

from pipeline_definition.types.type_registry import get_input_factory
from pipeline_definition.types.type_registry import get_step_factory

from pipeline_definition.types.input_type import InputFactory

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

    def translateYamlDoc(self, yamlDoc):

        #We need to know the inputs and outputs

        inputs = yamlDoc.get("inputs")
        steps = yamlDoc.get("steps")
        outputs = yamlDoc.get("outputs")

        if inputs is None:
            raise ValueError("No input?")

        if outputs is None:
            #raise ValueError("No output?")
            pass

        self.translateInputs(inputs)
        self.translateSteps(steps)

        return None

    def translateInputs(self, inputs ):
        for key, value in inputs.items():
            print("INPUT: " + key)

            inpFactory = get_input_factory( key )
            if ( inpFactory is None ):
                raise ValueError("No factory registered for input: " + key )

            val = inpFactory.emit()

            print(key, "=>", val)

        return None





    def translateSteps(self, steps ):

        for step in steps:
            #print("STEP: ", step, type(step))

            stepType = None
            if ( isinstance( step, dict) ):
                stepType = next( iter(step.keys()) )
            elif ( isinstance( step, str) ):
                stepType = step

            print("STEP: ", stepType)

            stepFactory = get_step_factory(stepType)

            if ( stepFactory is None ):
                raise ValueError("No factory registered for step: " + stepType )

            val = stepFactory.emit()

            print(stepType, "=>", val)


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

        tDoc = self.translateYamlDoc( doc )






