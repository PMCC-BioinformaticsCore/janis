class CwlIdLabelDoc:
    kID = "id"
    kLABEL = "label"
    kDOC = "doc"

class CwlHintsRequirements:
    kHINTS = "hints"
    kREQUIREMENTS = "requirements"


class Cwl(CwlIdLabelDoc):
    kCUR_VERSION = "v1.0"
    kCWL_VERSION = "cwlVersion"

    kCLASS = "class"

    class Class:
        kWORKFLOW = "Workflow"
        kCOMMANDLINETOOL = "CommandLineTool"
        kEXPRESSIONTOOL = "ExpressionTool"

    class Workflow(CwlIdLabelDoc, CwlHintsRequirements):

        kINPUTS = "inputs"
        kOUTPUTS = "outputs"
        kSTEPS = "steps"

        class Input(CwlIdLabelDoc):
            kTYPE = "type"
            kFORMAT = "format"
            kSTREAMABLE = "streamable"
            kDEFAULT = "default"
            kSECONDARY_FILES = "secondaryFiles"
            kINPUT_BINDING = "inputBinding"

            class InputBinding:
                kLOAD_CONTENTS = "loadContents"
                kPOSITION = "position"
                kPREFIX = "prefix"
                kSEPARATE = "separate"
                kITEM_SEPARATOR = "itemSeparator"
                kVALUE_FROM = "valueFrom"
                kSHELL_QUOTE = "shellQuote"

        class Output(CwlIdLabelDoc):
            kTYPE = "type"
            kSTREAMABLE = "streamable"
            kFORMAT = "format"
            kLINK_MERGE = "linkMerge"

            kSECONDARY_FILES = "secondaryFiles"
            kOUTPUT_SOURCE = "outputSource"
            kOUTPUT_BINDING = "outputBinding"

            class OutputBinding:
                kGLOB = "glob"
                kLOAD_CONTENTS = "loadContents"
                kOUTPUT_EVAL = "outputEval"

        class Step(CwlIdLabelDoc, CwlHintsRequirements):
            kRUN = "run"
            kSCATTER = "scatter"
            kSCATTER_METHOD = "scatterMethod"         # not supported
            kIN = "in"
            kOUT = "out"

            class StepInput:
                kID = "id"
                kSOURCE = "source"
                kLINK_MERGE = "linkMerge"
                kDEFAULT = "default"
                kVALUE_FROM = "valueFrom"

            class ScatterMethod:
                kDOT_PRODUCT = "dotproduct"
                kNESTED_CROSS_PRODUCT = "nested_crossproduct"
                kFLAT_CROSS_PRODUCT = "flat_crossproduct"

    class CommandLineTool(CwlIdLabelDoc, CwlHintsRequirements):
        kBASE_COMMAND = "baseCommand"

        kSTDIN = "stdin"
        kSTDERR = "stderr"
        kSTDOUT = "stdout"

        kSUCCESS_CODES = "successCodes"
        kTEMPORARY_FAIL_CODES = "temporaryFailCodes"
        kPERMANENT_FAIL_CODES = "permanentFailCodes"

        kINPUTS = "inputs"
        kOUTPUTS = "outputs"
        kARGUMENTS = "arguments"

        class Inputs(CwlIdLabelDoc):
            kTYPE = "type"
            kFORMAT = "format"
            kSTREAMABLE = "streamable"
            kDEFAULT = "default"
            kSECONDARY_FILES = "secondaryFiles"
            kINPUT_BINDING = "inputBinding"

            class InputBinding:
                kLOAD_CONTENTS = "loadContents"
                kPOSITION = "position"
                kPREFIX = "prefix"
                kSEPARATE = "separate"
                kITEM_SEPARATOR = "itemSeparator"
                kVALUE_FROM = "valueFrom"
                kSHELL_QUOTE = "shellQuote"

        class Outputs(CwlIdLabelDoc, CwlHintsRequirements):
            kTYPE = "type"
            kFORMAT = "format"
            kSTREAMABLE = "streamable"
            kDEFAULT = "default"
            kSECONDARY_FILES = "secondaryFiles"
            kOUTPUT_SOURCE = "outputSource"
            kOUTPUT_BINDING = "outputBinding"

            class OutputBinding(CwlIdLabelDoc):
                kGLOB = "glob"
                kLOAD_CONTENTS = "loadContents"
                kOUTPUT_EVAL = "outputEval"

        class Arguments:
            kLOAD_CONTENTS = "loadContents"
            kPOSITION = "position"
            kPREFIX = "prefix"
            kSEPARATE = "separate"
            kITEM_SEPARATOR = "itemSeparator"
            kVALUE_FROM = "valueFrom"
            kSHELL_QUOTE = "shellQuote"

    class Requirements:
        kCLASS = "class"
        kJAVASCRIPT = "InlineJavascriptRequirement"
        kSUBWORKFLOW = "SubworkflowFeatureRequirement"
        kSCATTER = "ScatterFeatureRequirement"
        kMULTIPLEINPUT = "MultipleInputFeatureRequirement"
        kSTEPINPUTEXPRESSION = "StepInputExpressionRequirement"
        kSCHEMADEF = "SchemaDefRequirement"
        kSOFTWARE = "SoftwareRequirement"
        kINITIALWORKDIR = "InitialWorkDirRequirement"

    class Primitives:
        kNULL = "null"
        kBOOLEAN = "boolean"
        kINT = "int"
        kLONG = "long"
        kFLOAT = "float"
        kDOUBLE = "double"
        kSTRING = "string"
        kFILE = "File"
        kDIRECTORY = "Directory"
        kARRAY = "array"
