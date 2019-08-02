"""
    This file declares the 'entry_points' that Janis will look for in its operations, and
    that your tool can bind to in order to provide additional functionality to Janis.

    > Check this awesome tutorial about Snek Inc for more information about entry-points:
        https://amir.rachum.com/blog/2017/07/28/python-entry-points/

    This is how all of janis works, the unix and bioinformatics tools are registered under
    all three of 'janis.extension', 'janis.types' and 'janis.tools'. janis-runner is just
    bound under 'janis.extension' to provide 'janis.runner' and CLI functionality.
"""


"""
    All extensions registered under this entrypoint will be added as an export to 
    janis.{yourextension}. Additionally, they will be added to the python sys.path 
    so you can just call'import janis.{yourextension}'.
"""
EXTENSIONS = "janis.extension"


"""
    All extensions that want to declare types should register a 'janis.types' extension.
    The shed (`janis.shed`) will look up this entry point and traverse any child modules
    looking for data types. 
    
    Note: Users can still access your type through a direct Python import (preferred)

"""
DATATYPES = "janis.types"

"""
    All extensions that want to provide tools should register a 'janis.tools' extension.
    The shed (`janis.shed`) will look up this entry point and traverse any child modules
    looking for tools. 
    
    Note: Users can still access your tool through a direct Python import (preferred)
"""
TOOLS = "janis.tools"
