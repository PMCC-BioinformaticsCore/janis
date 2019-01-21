# Documentation: https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md#types
import logging
from typing import Optional

logging.basicConfig(level=logging.INFO)
_LOGGER = logging.getLogger(__name__)


class PrimitiveType:

    kBoolean = "Boolean"
    kInt = "Int"
    kFloat = "Float"
    kString = "String"
    kFile = "File"

    DEF_TYPE = kString

    types = [kBoolean, kInt, kFloat, kFile, kString]

    def __init__(self, prim_type):
        if prim_type not in self.types:
            raise Exception("'{t}' was not a primitive type, expected one of {types}"
                            .format(t=prim_type, types=", ".join(self.types)))
        self._type = prim_type

    def get_string(self):
        return self._type

    @staticmethod
    def parse(prim_type):
        if prim_type not in PrimitiveType.types:
            return None
        return PrimitiveType(prim_type)


class ArrayType:

    kArray = "Array"

    def __init__(self, subtype, requires_multiple):

        self._subtype: WdlType = WdlType.parse_type(subtype, requires_type=True)
        self._requires_multiple: bool = requires_multiple

    def get_string(self):

        f = ArrayType.kArray + "[{t}]{quantifier}"

        if isinstance(self._subtype, list):
            return [f.format(t=t.get_string(), quantifier=("+" if self._requires_multiple else "")) for t in self._subtype]

        wd = self._subtype.get_string()
        if isinstance(wd, list) and len(wd) > 1:
            raise Exception("Internal error: unable to support array type with multiple subvalues")
        return f.format(
            t=wd[0] if isinstance(wd, list) else wd,
            quantifier=("+" if self._requires_multiple else "")
        )

    @staticmethod
    def parse(t: str, requires_multiple: Optional[bool]=None):

        requires_multiple = (requires_multiple is True) or False
        if t.endswith("+"):
            requires_multiple = True
            t = t[:-1]

        if not (t.startswith(ArrayType.kArray + "[") and t.endswith("]")):
            return None

        return ArrayType(t[6:-1], requires_multiple)

# class CompoundTypes(ArrayType, PrimitiveType):


class WdlType:

    postfix_quantifiers = [
        # Documentation: https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md#optional-parameters--type-constraints
        "?",    # means that the value is optional. if the value is missing, it will evaluate to the empty string.
        "+"     # can only be applied to Array types, the array is required to have one or more values in it
    ]

    types = [*PrimitiveType.types, ArrayType.kArray]

    def __init__(self, type_obj, optional=False):
        if not (isinstance(type_obj, PrimitiveType)
                or isinstance(type_obj, ArrayType)
                or isinstance(type_obj, WdlType)):
            raise Exception("Must initialise WdlType with PrimitiveType, ArrayType or WdlType")

        self._type = type_obj
        self.optional = optional

    def get_string(self):

        wd = self._type.get_string()
        if isinstance(wd, list):
            return [t + ("?" if self.optional else "") for t in wd]
        else:
            return wd + ("?" if self.optional else "")

    @staticmethod
    def parse_type(t, requires_type=True):
        """
        Will parse the type
        :param t:
        :param requires_type:
        :return:
        """
        t_orig = str(t)

        if not t:
            raise Exception("Must pass a value to parse_type")

        if isinstance(t, list):
            return [WdlType.parse_type(tt, requires_type) for tt in t]

        if isinstance(t, WdlType):
            return t
        if isinstance(t, PrimitiveType) or isinstance(t, ArrayType):
            return WdlType(t)

        if isinstance(t, str):
            optional_quantifier, multi_quantifier, t = WdlType.check_quantifiers(t)

            parse_attempt1 = PrimitiveType.parse(t)
            if parse_attempt1:
                if multi_quantifier: _LOGGER.warning("Ignoring extraneous multi_quantifier (+) on '{tp}'".format(tp=t_orig))
                return WdlType(parse_attempt1, optional=optional_quantifier)

            parse_attempt2 = ArrayType.parse(t, multi_quantifier)
            if parse_attempt2:
                return WdlType(parse_attempt2, optional=optional_quantifier)

        if requires_type:
            raise Exception("Couldn't pass '{t}'".format(t=t_orig))

        _LOGGER.warning(f"Returning None type for '{t_orig}'")
        return None

    @staticmethod
    def check_quantifiers(t: str):

        optional_quantifier = False
        multi_quantifier = False

        if t.endswith("?"):
            t = t[:-1]
            optional_quantifier = True

        if t.endswith("+"):
            t = t[:-1]
            multi_quantifier = True

        if not optional_quantifier and t.endswith("?"):
            t = t[:-1]
            optional_quantifier = True

        return optional_quantifier, multi_quantifier, t


Boolean = WdlType(PrimitiveType(PrimitiveType.kBoolean))
Int = WdlType(PrimitiveType(PrimitiveType.kInt))
Float = WdlType(PrimitiveType(PrimitiveType.kFloat))
File = WdlType(PrimitiveType(PrimitiveType.kFile))
String = WdlType(PrimitiveType(PrimitiveType.kString))
