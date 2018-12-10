class Cwl:
    kCUR_VERSION = "v1.0"

    kCWL_VERSION = "cwlVersion"

    kWORKFLOW = "Workflow"
    kTOOL = "CommandLineTool"

    class WORKFLOW:
        kCLASS = "class"
        kINPUTS = "inputs"
        kOUTPUTS = "output"
        kSTEPS = "steps"
        kHINTS = "hints"
        kID = "id"

        class STEP:
            kIN = "in"
            kOUT = "out"
            kRUN = "run"
            kSCATTER = "scatter"

    class COMMANDLINETOOL:
        kCLASS = "class"

    kREQUIREMENTS = "requirements"

    class REQUIREMENTS:
        kSCATTER_FEATURE = "ScatterFeatureRequirement"
        kSUBWORKFLOW_FEATURE = "SubworkflowFeatureRequirement"
        kJAVASCRIPT = "InlineJavascriptRequirement"

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

