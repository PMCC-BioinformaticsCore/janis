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
from pipeline_definition.types.step_type import Step
from pipeline_definition.types.InputStep import InputStep

class PipelineTranslatorException(Exception):
    pass

class PipelineTranslator:
    def __init__(self):
        self.__inputStep = None

    def __dumpYaml(self, doc ):
        # Diagnostic - what have we got?
        print("YAML DOC [")
        print(yaml.dump(doc))
        print("] END YAML DOC")

    def __dumpSchema(self, schema ):
        print("PDX SCHEMA [")
        #print(schema)
        print(json.dumps(schema, indent=4))
        print("] END PDX SCHEMA")


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
                stepType = Step.selectTypeNameFrom( meta )

            print("Processing STEP: ",id, " - ", stepType)

            stepFactory = get_step_factory( stepType )
            if ( stepFactory is None ):
                raise ValueError("No factory registered for step: " + stepType )

            stepObj = stepFactory.buildFrom( dict([ (id, meta) ]) )

            stepObj.validateInputOuputSpec()

            pipelineSteps.append( stepObj )


        return pipelineSteps


    def __createWorkflowGraph(self, pipelineSteps, workflowInputSet, workflowOutputSet ):
        workGraph = nx.MultiDiGraph()

        #Lets create the input step - the start step of every workflow that produces the workflow inputs as step output
        self.__inputStep = InputStep(workflowInputSet)
        pipelineSteps.insert(0, self.__inputStep)

        #Now lets put all the steps as node. Then we will establish the edges
        print("Graph Construction: adding nodes for steps.")
        for step in pipelineSteps:
            stepCtx = StepContext(step)
            workGraph.add_node(step, ctx=stepCtx)
            print("Added node for step", step.id(),"[",step.tag(),"]")

        print("Graph Construction: Processing tags.")
        #Now lets put the edges to indicate execution order as indicated by 'tag'
        #Convention is all steps belonging to a tag are exceuted in the seq they have been specified
        #so creating the concept of 'threads' / 'branches'

        #Create a dict of tags vocabulary that have been used in the workflow descriptin
        tagMap = dict()
        for step in pipelineSteps:
            tag = step.tag()
            if tag not in tagMap:
                tagMap[tag] = None

        print("TAG MAP:", str(tagMap))

        #Lets stitch the DAG based on tags assigned to each step

        for step in pipelineSteps:
            print("Processing step ", step.id(), "[", step.tag(), "]")

            #tag specification in step
            stag = step.tag()
            #print("Step specifies tag:", stag)

            if ( stag == self.__inputStep.tag() ):
                tagMap[stag] = self.__inputStep
                continue

            #Have we already seen a step in that thread/tag?
            lastNode = tagMap.get(stag)
            tagMap[stag] = step

            if lastNode is None:
                lastNode = self.__inputStep

            workGraph.add_edge(lastNode, step, type="branch", tag=stag)
            print("Added edge [",lastNode.id(),"]->[",step.id(),"] for tag", stag)

        #Next pass is about input dependency of a step on outputs form other branches and steps
        print("Graph Construction: Processing input dependency")
        for step in pipelineSteps:
            print("Processing step ", step.id(), "[", step.tag(),"]")

            if step == self.__inputStep:
                print("Input step never have dependency.")
                continue

            #Lets get the dependency list of the step
            depends = self.dependencyListOf( step )
            if not depends:
                print("Step has no dependency.")
                continue

            for dependency in depends:
                print("Step has dependecny:", dependency)
                self.addDependencyTo(step, dependency, workGraph)

        print("Graph Construction: Polulating context for steps.")
        #Now we have a graph that has edges for flow and dependency
        #Now we need to establish the context of each step
        ctxAttrMap = nx.get_node_attributes(workGraph, 'ctx')
        typeAttrMap = nx.get_edge_attributes(workGraph, 'type')
        self.__populateSetpContext(workGraph, self.__inputStep, None, ctxAttrMap, typeAttrMap )

        return workGraph

    def __populateSetpContext(self, workGraph, step, prevStep, ctxAttrMap, typeAttrMap ):
        print("Populating context for step:", step.id(),"in branch [", step.tag(),"]" )
        stepCtx = ctxAttrMap.get(step)
        if not stepCtx:
            raise RuntimeError("Missing step context in graph. Graph integrity fail.")

        if prevStep:
            prevCtx = ctxAttrMap.get(prevStep)
            if not prevCtx:
                raise RuntimeError("Missing step context in graph. Graph integrity fail.")

            if prevStep.tag() != InputStep.inputSteptagName() and step.tag() != prevStep.tag():
                raise RuntimeError("Branch tag mismatch during context population.")

            stepCtx.inheritContextOfBranch(prevCtx)

        edges = nx.edges(workGraph, step)
        if not edges:
            return

        print("Step [",step.id(),"] in branch [",step.tag(),"] has", len(edges), "edges")

        #First pass to process the input dependency
        for edge in edges:
            if edge[0] != step:
                raise RuntimeError("Trouble in edge finding")
            edgeType = str(typeAttrMap[(edge[0], edge[1], 0)])
            if edgeType != 'dependency':
                continue
            dependentOn = edge[1]
            dependentOnCtx = ctxAttrMap[dependentOn]
            print("Step",step.id()+"["+step.tag()+"] has dependency on step",dependentOn.id(),"[",dependentOn.tag(),"]")
            stepCtx.addDependencyContextFrom(dependentOnCtx)

        stepCtx.print()

        #Sencod pass to process the flow
        for edge in edges:
            if edge[0] != step:
                raise RuntimeError("Trouble in edge finding")
            edgeType = str(typeAttrMap[(edge[0], edge[1], 0)])
            if edgeType != 'branch':
                continue
            self.__populateSetpContext( workGraph, edge[1], step, ctxAttrMap, typeAttrMap)

        return



    def addDependencyTo(self, step, dependencySpec, workGraph):

        if not step:
            return

        if not dependencySpec:
            return

        if not workGraph:
            return

        targetTag = dependencySpec.get('tag')
        targetStep = dependencySpec.get('step')

        #Lets start on the graph by tag
        sourceStep = self.__inputStep
        edgeMapforTags = nx.get_edge_attributes(workGraph, 'tag')

        lastStep = self.lastStepInTag(targetTag, sourceStep, workGraph, edgeMapforTags)
        if lastStep:
            print(step.id(),"in branch [", step.tag(), "] has input dependency on step", lastStep.id(), "in branch [",targetTag,"]")
            workGraph.add_edge(step, lastStep, type="dependency")

        return

    def lastStepInTag(self, tag, root, workGraph, edgeMapforTags ):
        sourceStep = root
        edges = nx.edges(workGraph, sourceStep)
        if not edges:
            #there are no edges so this is the last step
            return sourceStep

        for edge in edges:
            edgeTag = str(edgeMapforTags[(edge[0], edge[1], 0)])
            if edgeTag == tag:
                return self.lastStepInTag( tag, edge[1], workGraph, edgeMapforTags)
        return sourceStep


    def dependencyListOf(self, step ):
        if not step:
            return None

        stepRequires = step.requires()
        if not stepRequires:
            return None

        dependencyList = None
        for requirement in stepRequires:
            requirementName = requirement[Step.STR_ID]
            #print("Process STEP REQUIREMENT:", requirementName)
            requirementValue = step.providedValueForRequirement(requirementName)
            #print("Input value:", requirementValue)

            if not requirementValue:
                continue

            dependencySpec = self.dependencySpecFrom( requirementValue )

            if not dependencyList:
                dependencyList = list()
            dependencyList.append(dependencySpec)

        return dependencyList


    def dependencySpecFrom(self, requirementValue ):

        tag = None
        step = None
        output = None

        parts = requirementValue.split(".")

        head = parts[0]
        if head.startswith("#"):
            tag = head[1:]

        return {
            'tag' : tag,
            'step' : step,
            'output' : output
        }


    def __dumpGraph(self, workGraph):
        #tree = json_graph.tree_data(workGraph, self.__root, attrs={'children': 'next', 'id': 'step'})
        tree = json_graph.node_link_data(workGraph,{'link': 'flow', 'source': 'step', 'target': 'target'})
        print("Workflow Graph: [")
        print(tree)
        #jsonDoc = json.dumps(tree, indent=4)
        #print(jsonDoc)

        print("] End Workflow Graph")


    def translatePipeline( self, pipelineSteps, globalInputSet, globalOutputSet ):

        workGraph = self.__createWorkflowGraph( pipelineSteps, globalInputSet, globalOutputSet )
        self.__dumpGraph( workGraph )

        jsonDoc = self.translateWorkflowToJson(workGraph)

        prettyJsonText = json.dumps(jsonDoc, indent=4)

        return prettyJsonText

    def translateWorkflowToJson(self, workGraph):
        print("Generating JSON description for workflow")
        ctxAttrMap = nx.get_node_attributes(workGraph, 'ctx')
        typeAttrMap = nx.get_edge_attributes(workGraph, 'type')
        workflow = self.describeWorkflow(workGraph, self.__inputStep, typeAttrMap, ctxAttrMap)
        jsonDoc = {}
        jsonDoc["workflow"] = workflow
        print("Done generating JSON description for workflow")
        return jsonDoc

    def describeWorkflow(self, workGraph, inputStep, typeAttrMap, ctxAttrMap ):

        doc = {
        }

        stepCtx = ctxAttrMap[inputStep]
        provides = stepCtx.provides()

        doc['inputs'] = inputStep.provides()


        flow = {}
        edges = nx.edges(workGraph, inputStep)
        if not edges:
            return

        for edge in edges:
            if edge[0] != inputStep:
                raise RuntimeError("Trouble in edge finding")
            edgeType = str(typeAttrMap[(edge[0], edge[1], 0)])
            if edgeType == 'branch':
                branchTag = edge[1].tag()
                branchDesc = self.__branchDescriptionFrom(edge[1])
                flow[ branchTag ] = branchDesc



        doc['flow'] = flow

        return doc

    def __branchDescriptionFrom( self, step ):



        return {}


    def translateYamlDoc(self, yamlDoc):
        #Create a in memory instances for all the inputs, steps and outputs

        inputs = yamlDoc.get("inputs")
        steps = yamlDoc.get("steps")
        outputs = yamlDoc.get("outputs")

        if inputs is None:
            raise ValueError("No input?")

        #if outputs is None:
        #    raise ValueError("No output?")

        workflowInputSet = self.buildInputs(inputs)
        workflowOutputSet = self.buildOutputs(outputs)
        pipelineSteps = self.buildSteps(steps )

        #Now translate the workflow steps
        doc = self.translatePipeline( pipelineSteps, workflowInputSet, workflowOutputSet )

        return doc

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
            outputDoc = self.translateYamlDoc( doc )
            print("Translation Output: [")
            print(outputDoc)
            print("]")

            # Save it
            #with open( outfilePath, "w" ) as ofile:
            #    ofile.write( tDoc )

