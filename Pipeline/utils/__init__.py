from typing import Dict


def first_value(d: Dict):
    return next(iter(d.values()))


def convert_expression_to_wdl(expression):
    import re

    # Some examples
    # $(inputs.filename) -> "{filename}"
    # $(inputs.filename).out -> "{filename}.out"
    # randomtext.filename -> "randomtext.filename

    if not expression:
        return expression

    if not isinstance(expression, str):
        raise Exception(f"Expected expression of type string, got: {expression} ({type(expression)})")

    r1 = re.compile(r"\$\(inputs\..*\)")
    m1 = r1.match(expression)
    if m1:
        return f'"{{{m1.group()[9:-1]}}}{expression[m1.span()[1]:]}"'

    r2 = re.compile(r"\$\(\..*\)")
    m2 = r2.match(expression)
    if m2:
        return f'"{{{m2.group()[9:-1]}}}{expression[m1.span()[1]:]}"'

    return f'"{expression}"'
