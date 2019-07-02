import warnings


class PipelineTranslatorException(Exception):
    pass


class DuplicateLabelIdentifier(Exception):
    pass


class InvalidNodeIdentifier(Exception):
    pass


class NodeNotFound(Exception):
    pass


class NotFoundException(Exception):
    pass


class InvalidInputsException(Exception):
    pass


class InvalidStepsException(Exception):
    pass


class TooManyArgsException(Exception):
    pass


class IncorrectArgsException(Exception):
    pass


class InvalidByProductException(Exception):
    pass


class ConflictingArgumentsException(Exception):
    pass


def deprecated(message):
    def deprecated_decorator(func):
        def deprecated_func(*args, **kwargs):
            warnings.warn(
                "{} is a deprecated function. {}".format(func.__name__, message),
                category=DeprecationWarning,
                stacklevel=2,
            )
            warnings.simplefilter("default", DeprecationWarning)
            return func(*args, **kwargs)

        return deprecated_func

    return deprecated_decorator
