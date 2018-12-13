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

    class CLASS:
        kWORKFLOW = "Workflow"
        kCOMMANDLINETOOL = "CommandLineTool"
        kEXPRESSIONTOOL = "ExpressionTool"

    class WORKFLOW(CwlIdLabelDoc, CwlHintsRequirements):

        kINPUTS = "inputs"
        kOUTPUTS = "outputs"
        kSTEPS = "steps"

        class INPUT(CwlIdLabelDoc):
            kTYPE = "type"
            kFORMAT = "format"
            kSTREAMABLE = "streamable"
            kDEFAULT = "default"
            kSECONDARY_FILES = "secondaryFiles"
            kINPUT_BINDING = "inputBinding"

            class INPUT_BINDING:
                kLOAD_CONTENTS = "loadContents"
                kPOSITION = "position"
                kPREFIX = "prefix"
                kSEPARATE = "separate"
                kITEM_SEPARATOR = "itemSeparator"
                kVALUE_FROM = "valueFrom"
                kSHELL_QUOTE = "shellQuote"

        class OUTPUT(CwlIdLabelDoc):
            kTYPE = "type"
            kSTREAMABLE = "streamable"
            kFORMAT = "format"
            kLINK_MERGE = "linkMerge"

            kSECONDARY_FILES = "secondaryFiles"
            kOUTPUT_SOURCE = "outputSource"
            kOUTPUT_BINDING = "outputBinding"

            class OUTPUT_BINDING:
                kGLOB = "glob"
                kLOAD_CONTENTS = "loadContents"
                kOUTPUT_EVAL = "outputEval"



        class STEP(CwlIdLabelDoc, CwlHintsRequirements):
            kRUN = "run"
            kSCATTER = "scatter"
            kSCATTER_METHOD = "scatterMethod"         # not supported
            kIN = "in"
            kOUT = "out"

            class STEP_INPUT:
                kID = "id"
                kSOURCE = "source"
                kLINK_MERGE = "linkMerge"
                kDEFAULT = "default"
                kVALUE_FROM = "valueFrom"

            class SCATTER_METHOD:
                kDOT_PRODUCT = "dotproduct"
                kNESTED_CROSS_PRODUCT = "nested_crossproduct"
                kFLAT_CROSS_PRODUCT = "flat_crossproduct"

    class COMMANDLINETOOL(CwlIdLabelDoc, CwlHintsRequirements):
        pass

    class REQUIREMENTS:
        kCLASS = "class"
        kJAVASCRIPT = "InlineJavascriptRequirement"
        kSUBWORKFLOW = "SubworkflowFeatureRequirement"
        kSCATTER = "ScatterFeatureRequirement"
        kMULTIPLEINPUT = "MultipleInputFeatureRequirement"
        kSTEPINPUTEXPRESSION = "StepInputExpressionRequirement"
        kSCHEMADEF = "SchemaDefRequirement"
        kSOFTWARE = "SoftwareRequirement"
        kINITIALWORKDIR = "InitialWorkDirRequirement"

    class PRIMITIVES:
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

