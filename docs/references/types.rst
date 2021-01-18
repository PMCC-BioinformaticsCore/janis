Types
#############

Janis has a typing system that allows pipeline developers to be confident when they connect components together. In Janis, there are 6 fundamental types, all of which are inheritable. These fundamental types require an equivalent in all translation targets (eg: ``janis.String`` must be to a string CWL and WDL type).

- :class:`janis.File`
- :class:`janis.String`
- :class:`janis.Int`
- :class:`janis.Float`
- :class:`janis.Boolean`
- :class:`janis.Directory`

To provide greater context to your tools and workflows, you can create a custom type that inherits from all of these fundamental types.

Often you can just use the uninstantiated type for annotation. All types can be made optional by instantiating the type wth ``(optional=True)``. By default types are required (non-optional) and inputs that have a default are made optional.

Patterns
**********

The janis-patterns repository has an example of how python types relate to Janis types: https://github.com/PMCC-BioinformaticsCore/janis-patterns/blob/master/types/pythontypes.py

File
********

A file in Janis is an object, and not a path. By annotating your type as a file, the execution engine will make the referenced file available within your execution environment. You should simply consider that the File type is a reference to some arbitrary place on disk and not in a predictable location.

Most often you will not directly instantiate a File as you'll use some domain specific file type (such as ``TarFile`` or `Bam``).

.. autoclass:: janis.File

If you do not correctly annotate your type to be a :class:`janis.File` (with appropriate index files), it may not be localised into your execution environment.

Secondary / Accessory files
============================

Janis inherits the concept of secondary files from the Common Workflow Language. When creating a subclass of a File, you should include the

.. code-block:: python

   @staticmethod
   def secondary_files() -> Optional[List[str]]:
       return None

And return a list of strings with the following syntax:

- If string begins with one or more caret ``^`` characters, for each caret, remove the last file extension from the path (the last period ``.`` and all following characters). If there are no file extensions, the path is unchanged.
- Append the remainder of the string to the end of the file path.

.. note::

   It's important that if you're directly instantiating the File class, or subclassing and calling ``super().__init__`` that you MUST include the expected file extension, otherwise secondary files will be not correctly collected in WDL.


Typing system
==============

The typing system allows for you to pass an inherited type to one of its super classes. This means that if your type *inherits* from ``File``, you can safely pass your inherited type. When passing inherited file types with secondary files, Janis does it's best to take the subset of files from the receiving type. It's up to the developer's of new data types to ensure that all secondary files are included from each subtype.

String
********

A ``String`` is a type in Janis. When annotating a type, you can also use the inbuilt Python ``str``, which will be turned into a required (non-optional) :class:`janis.String`.

.. autoclass:: janis.String

The :class:`janis.Filename` class is a special type that inherits from string which allows for automatic filename generation.

Int
********

An ``Int`` is a type in Janis. When annotating a type, you can also use the inbuilt Python ``int``, which will be turned into a required (non-optional) :class:`janis.Int`.

.. autoclass:: janis.Int

The range of an integer is the minimum of the specification you use, here are some values:

========  =========  ===============
Language  Min        Max
========  =========  ===============
CWL       ``-2^32``  ``2^32``
WDL       ``-2^63``  ``2^63``
Python    →          ``sys.maxsize``
JSON      ``*``      ``*``
YAML      ``*``      ``*``
========  =========  ===============

``*`` - JSON and YAML max ranges may depend on your parser.

Float
********

A ``Float`` is a type in Janis. When annotating a type, you can also use the inbuilt Python ``float``, which will be turned into a required (non-optional) :class:`janis.Float`.

.. autoclass:: janis.Float

The range of an integer is the minimum of the specification you use, here are some values:

========  =========  =======================
Language  Min        Max
========  =========  =======================
CWL       →          finite 32-bit IEEE-754
WDL       →          finite 64-bit IEEE-754
Python    →          ``sys.float_info``
JSON      ``*``      ``*``
YAML      ``*``      ``*``
========  =========  =======================

``*`` - JSON and YAML max ranges may depend on your parser.

Boolean
********

A ``Boolean`` is a type in Janis. When annotating a type, you can also use the inbuilt Python ``bool``, which will be turned into a required (non-optional) :class:`janis.Boolean`.

.. autoclass:: janis.Boolean

Filename
********

A ``Filename`` is an inherited type in Janis with the purpose to stand in for filename generation. By default it is always optional, and there are a number of ways you can customise it.

.. autoclass:: janis.Filename

- Prefix (defaults to "generated")
- Suffix
- Guid (defaults to evaluating ``str(uuid.uuid1())`` Overridable if you want a similar structure to your output files + some extension.
- Extension (this is very important if collecting output files), you must include the period (``.``).

Format: ``$prefix-$guid-$suffix.bam``

The intention was for each run and scatter to have a different generated value for each instance of generated filename. However, at the moment technical difficulties has delayed this implementation.

