import re


class Validators:
    identifier_regex = r"[a-zA-Z][a-zA-Z0-9_]+\Z"
    compiled_identifier = re.compile(identifier_regex)

    @staticmethod
    def validate_identifier(identifier) -> bool:
        a = Validators.compiled_identifier.match(identifier)
        return a is not None
