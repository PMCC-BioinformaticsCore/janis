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

        inputSet = self.translateInputs(inputs)
        outputSet = self.translateOutputs(outputs)
        self.translateSteps(steps, inputSet, outputSet )

        return "TX Output"

    def translateInputs(self, inputs ):
        for key, value in inputs.items():
            print("INPUT: " + key)

            inpFactory = get_input_factory( key )
            if ( inpFactory is None ):
                raise ValueError("No factory registered for input: " + key )

            val = inpFactory.emit()
            print(key, "=>", val)

            inputObj = inpFactory.build( dict([ (key,value) ]) )
        return None

    def translateOutputs(self, outputs ):
        return None

    def translateSteps(self, steps, inputSet, outputSet ):

        for step in steps:
            #print("STEP: ", step, type(step))

            key = None
            value = None
            if ( isinstance( step, dict) ):
                key = next( iter(step.keys()) )
                value = next( iter(step.values()) )
            elif ( isinstance( step, str) ):
                key = step
                value = None

            print("STEP: ", key)

            stepFactory = get_step_factory(key)

            if ( stepFactory is None ):
                raise ValueError("No factory registered for step: " + key )

            val = stepFactory.emit()
            print(key, "=>", val)

            stepObj = stepFactory.build( dict([ (key,value) ]) )



    def translate(self, pdfile, outfile=None, overwriteOutfile=False ):
        pdfilePath = os.path.abspath(pdfile)
        print("Using PD file: " + pdfilePath)

        if not os.path.isfile(pdfilePath):
            raise ValueError("PD file does not exists.")

        if outfile is None:
            outfile = pdfile + ".pdx"
        outfilePath = os.path.abspath(outfile)

        print("Using Output file: " + outfilePath)
        if overwriteOutfile is False:
            if os.path.isfile(outfilePath) :
                raise ValueError("Outfile already exists.")

            if os.path.isdir(outfilePath) :
                raise ValueError("Directoiry already exists with the name of outfile.")

        with open(pdfile, 'r') as ifile:
            # Validate YAML syntax
            doc = yaml.load(ifile)

            # Do schema validation
            self.validateSchema( doc )

            tDoc = self.translateYamlDoc( doc )
            print("Output: ", tDoc)

            with open( outfilePath, "w" ) as ofile:
                ofile.write( tDoc )

