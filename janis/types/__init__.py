from abc import ABC

from janis.utils.bracketmatching import get_keywords_between_braces
from janis.utils.errors import (
    TooManyArgsException,
    IncorrectArgsException,
    InvalidByProductException,
    ConflictingArgumentsException,
)
from janis.utils.logger import Logger


class Selector(ABC):
    pass


class InputSelector(Selector):
    def __init__(self, input_to_select, use_basename=None):
        # maybe worth validating the input_to_select identifier
        self.input_to_select = input_to_select
        self.use_basename = use_basename

    def to_string_formatter(self):
        kwarg = {self.input_to_select: self}
        return StringFormatter(f"{{{self.input_to_select}}}", **kwarg)

    def __add__(self, other):
        return self.to_string_formatter() + other


class WildcardSelector(Selector):
    def __init__(self, wildcard, first_element=None):
        self.wildcard = wildcard
        self.first_element = first_element


class MemorySelector(InputSelector):
    def __init__(self):
        super().__init__("runtime_memory")


class CpuSelector(InputSelector):
    def __init__(self, default=1):
        super().__init__("runtime_cpu")
        self.default = default


class StringFormatter(Selector):
    def __init__(self, format, **kwargs):
        self._format = format

        keywords, balance = get_keywords_between_braces(self._format)

        if balance > 0:
            Logger.warn(
                "There was an imbalance of braces in the string _format, this might cause issues with concatenation"
            )

        skwargs = set(kwargs.keys())

        if not keywords == skwargs:
            # what's the differences
            if not keywords.issubset(skwargs):
                raise IncorrectArgsException(
                    "The _format required additional arguments to be provided by "
                    "**kwargs, requires the keys:" + ", ".join(keywords - skwargs)
                )
            else:
                raise TooManyArgsException(
                    "The **kwargs contained unrecognised keys: "
                    + ", ".join(skwargs - keywords)
                )

        self.kwargs = kwargs

    resolved_types = [str, int, float]

    def resolve_with_resolved_values(self, **resolved_values):

        s1 = set(self.kwargs.keys())
        actual_keys, _ = get_keywords_between_braces(self._format)
        if s1 != actual_keys:
            diff = (actual_keys - s1).union(s1 - actual_keys)

            raise Exception(
                "The format for the string builder has changed since runtime, or an internal error has"
                " occurred. The following keys did not appear in both sets: "
                + ", ".join(diff)
            )

        s2 = set(resolved_values.keys())

        missing_keys = s1 - s2
        if len(missing_keys) > 0:
            raise IncorrectArgsException(
                "There were missing parameters when formatting string: "
                + ", ".join(missing_keys)
            )

        unresolved_values = [
            f"{r} ({type(resolved_values[r]).__name__})"
            for r in resolved_values
            if not any(
                isinstance(resolved_values[r], t)
                for t in StringFormatter.resolved_types
            )
        ]
        if len(unresolved_values) > 0:
            raise ValueError(
                "There were unresolved parameters when formatting string: "
                + ", ".join(unresolved_values)
            )

        retval = self._format
        for k in resolved_values:
            retval = retval.replace(f"{{{k}}}", str(resolved_values[k]))
        return retval

    def __radd__(self, other):
        return StringFormatter(other) + self

    def __add__(self, other):
        if isinstance(other, str):
            # check if it has args in it
            keywords = get_keywords_between_braces(other)
            if len(keywords) > 0:
                invalidkwargs = [k for k in self.kwargs if k not in self.kwargs]
                if len(invalidkwargs) > 0:
                    raise InvalidByProductException(
                        f"The string to be concatenated contained placeholder(s) ({', '.join(invalidkwargs)})"
                        f"that were not in the original StringFormatter"
                    )
            return self._create_new_formatter_from_strings_and_args(
                [self._format, other], **self.kwargs
            )

        elif isinstance(other, InputSelector):
            return self + other.to_string_formatter()

        elif isinstance(other, StringFormatter):
            # check if args overlap and they're different
            s1 = set(self.kwargs.keys())
            s2 = set(other.kwargs.keys())
            intersection = s1.intersection(s2)

            if len(intersection) > 0:
                not_same_args = [
                    k for k in intersection if self.kwargs[k] != other.kwargs[k]
                ]
                if len(not_same_args) > 0:
                    raise ConflictingArgumentsException(
                        f"Couldn't concatenate formats as there keys ({', '.join(not_same_args)}) "
                        f"that were not equal between formatters "
                    )

            # yeah we sweet
            new_args = {**self.kwargs, **other.kwargs}
            return StringFormatter._create_new_formatter_from_strings_and_args(
                [self._format, other._format], **new_args
            )

    @staticmethod
    def _create_new_formatter_from_strings_and_args(strings: [str], **kwargs):
        new_format = "".join(strings)
        try:
            return StringFormatter(new_format, **kwargs)
        except IncorrectArgsException as e:
            new_params = set(
                get_keywords_between_braces(new_format)[0] - set(kwargs.keys())
            )
            raise InvalidByProductException(
                "Joining the input files (to '{new_format}') created the new params: "
                + ", ".join(new_params)
            )
