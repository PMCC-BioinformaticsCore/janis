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
from pprint import pprint

from pipeline_definition.types.schema import schema

from pipeline_definition.types.type_registry import get_input_factory
from pipeline_definition.types.type_registry import get_step_factory

from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.utils.StepContext import StepContext

import networkx as nx
from networkx.readwrite import json_graph

class PipelineTranslatorException(Exception):
    pass

class PipelineTranslator:
    def __init__(self):
        self.__root = "start"

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

        sch = schema()
        self.__dumpSchema( sch )

        v = Validator(sch)

        validationSuccess = v.validate(yamlDoc)

        if not validationSuccess:
            msg = "ERROR! Pipeline definition document validation failed."
            print( msg )
            if v.errors is not None:
                raise ValueError( msg, v.errors )


    def buildInputs(self, inputs ):

        inputSet = list()

        for id, meta in inputs.items():

            if (isinstance(meta, str)):
                inputType = meta
                meta = dict([(inputType, None)])
            elif (isinstance(meta, dict)):
                inputType = next(iter(meta.keys()))

            print("Processing INPUT: ", id, " - ", inputType)

            inpFactory = get_input_factory( inputType )
            if ( inpFactory is None ):
                raise ValueError("No factory registered for input: " + inputType )

            inputObj = inpFactory.buildFrom( dict([ (id,meta) ]) )

            inputSet.append(inputObj)

            print("\n")

        return inputSet

    def buildOutputs(self, outputs ):
        return None

    def buildSteps(self, steps):

        pipelineSteps = list()

        for step in steps:

            id = next( iter(step.keys()) )
            meta = next( iter(step.values()) )

            if (isinstance(meta, str)):
                stepType = meta
                meta = dict([(stepType, None)])
            elif (isinstance(meta, dict)):
                stepType = next(iter(meta.keys()))

            print("Processing STEP: ",id, " - ", stepType)

            stepFactory = get_step_factory( stepType )
            if ( stepFactory is None ):
                raise ValueError("No factory registered for step: " + stepType )

            stepObj = stepFactory.buildFrom( dict([ (id, meta) ]) )

            pipelineSteps.append( stepObj )

            print("\n")

        return pipelineSteps

    def __createWorkflowGraph(self, pipelineSteps ):

        # create a graph of steps in pipeline
        # How many parallel threads we have? It is dictated by "tag".
        # During building the graph we need to keep track of the last step seen in a thread

        tagRootMap = dict()
        tagMap = dict()
        for step in pipelineSteps:
            tag = step.tag()
            if tag not in tagMap:
                tagMap[tag] = None

        print("TAG MAP:", str(tagMap))

        workGraph = nx.MultiDiGraph()

        workGraph.add_node( self.__root )

        # Pass one is to stitch the steps
        for step in pipelineSteps:
            stag = step.tag()

            # Have we already seen a step in that thread/tag?
            lastNode = tagMap[stag]

            # current step needs to be recorded as last node seen in the tag/thread
            tagMap[stag] = step

            #We create a context object for each step that will be decorated and use in later
            stepCtx = StepContext(step)
            #workGraph.add_node(step, attr_dict={ 'ctx' : stepCtx })
            workGraph.add_node(step, ctx=stepCtx )

            if lastNode is None:
                workGraph.add_edge(self.__root, step)
            else:
                workGraph.add_edge(lastNode, step)

        return workGraph

    def __dumpGraph(self, workGraph):
        tree = json_graph.tree_data(workGraph, self.__root, attrs={'children': 'next', 'id': 'step'})
        print("Workflow Graph: [\n")
        print(tree)
        #jsonDoc = json.dumps(tree, indent=4)
        #print(jsonDoc)

        print("] End Workflow Graph\n")


    def __populateStepContext(self, workGraph, step, prevstep, globalInputSet, globalOutputSet, contextualInputSet ):
        print("Populating context of step:", step.id() )

        ctxAttrs = nx.get_node_attributes(workGraph, 'ctx')
        stepCtx = ctxAttrs[step]
        #print("Step CTX:", stepCtx)

        stepCtx.setPrevStep( prevstep )
        stepCtx.setGlobalInputSet(globalInputSet)
        stepCtx.setGlobalOutputSet(globalOutputSet)

        if prevstep is not None:
            #Output from previous step can be input for current step
            stepOutputs =  prevstep.provides()
            if stepOutputs is not None:
                contextualInputSet[prevstep.id()] = stepOutputs
                stepCtx.setContextualInputSet(contextualInputSet)

        stepCtx.print()

        edges = nx.edges(workGraph, step)
        for e in edges:
            #print("Processing EDGE:",e)
            nstep = e[1]
            self.__populateStepContext(workGraph, nstep, step, globalInputSet, globalOutputSet, contextualInputSet)


    def __populateContexts(self, workGraph, globalInputSet, globalOutputSet ):

        #steps = nx.all_neighbors(workGraph, self.__root)
        edges = nx.edges(workGraph, self.__root)
        contextualInputSet = {}

        for e in edges:
            #print("Processing EDGE:",e)
            step = e[1]
            self.__populateStepContext(workGraph, step, None, globalInputSet, globalOutputSet, contextualInputSet )



    def translatePipeline( self, pipelineSteps, globalInputSet, globalOutputSet ):

        workGraph = self.__createWorkflowGraph( pipelineSteps )
        self.__dumpGraph( workGraph )

        self.__populateContexts(workGraph, globalInputSet, globalOutputSet )

        return None

    def translateYamlDoc(self, yamlDoc):
        #Create a in memory instances for all the inputs, steps and outputs

        inputs = yamlDoc.get("inputs")
        steps = yamlDoc.get("steps")
        outputs = yamlDoc.get("outputs")

        if inputs is None:
            raise ValueError("No input?")

        if outputs is None:
            #raise ValueError("No output?")
            pass

        globalInputSet = self.buildInputs(inputs)
        globalOutputSet = self.buildOutputs(outputs)
        pipelineSteps = self.buildSteps(steps )

        #Now translate the workflow steps
        txDoc = self.translatePipeline( pipelineSteps, globalInputSet, globalOutputSet )

        return txDoc

    def translate(self, pdfile, outfile=None, overwriteOutfile=False ):
        pdfilePath = os.path.abspath(pdfile)
        print("Using PD file: " + pdfilePath)

        #Check if the file to translate exists
        if not os.path.isfile(pdfilePath):
            raise ValueError("Specified PD file does not exist.")

        #The output file name, if not explicitly specified, a convention is used
        if outfile is None:
            outfile = pdfile + ".pdx"
        outfilePath = os.path.abspath(outfile)
        print("Using Output file: " + outfilePath)

        #If the output file path already exists as directory?
        if os.path.isdir(outfilePath):
            raise ValueError("Directoiry already exists with the name of outfile.")

        #Otherwise can we overwrite?
        if overwriteOutfile is False:
            if os.path.isfile(outfilePath) :
                raise ValueError("Outfile already exists.")

        #All good so lets starts the translation
        with open(pdfile, 'r') as ifile:
            # Validate YAML syntax
            doc = yaml.load(ifile)

            # Do schema validation
            self.validateSchema( doc )

            #Doc is OK, lets translate
            tDoc = self.translateYamlDoc( doc )
            print("Output: ", tDoc)

            # Save it
            #with open( outfilePath, "w" ) as ofile:
            #    ofile.write( tDoc )

