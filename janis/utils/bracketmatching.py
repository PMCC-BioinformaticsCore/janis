from keyword import iskeyword
from typing import Tuple

from janis.utils.logger import Logger


# https://stackoverflow.com/a/36331242/2860731
def variable_name_validator(x: str):
    return x.isidentifier() and not iskeyword(x)


def get_keywords_between_braces(
    text, argument_validator=variable_name_validator
) -> Tuple[set, int]:
    counter = 0
    highest_level = -1
    start_idx = None
    matches = set()

    for i in range(len(text)):

        char = text[i]
        if char == "{":
            counter += 1
            highest_level = max(highest_level, counter)
            if start_idx is None:
                start_idx = i
        elif char == "}" and counter > 0:
            counter -= 1

            if start_idx is not None and counter == 0:
                match = text[start_idx + 1 : i]
                if highest_level > 1:
                    Logger.log("Skipping match: " + match)
                elif argument_validator is not None and not argument_validator(match):
                    Logger.log("Match was rejected by validator: " + match)
                else:
                    Logger.log("Recognised match: " + match)
                    matches.add(match)
                highest_level = -1
                start_idx = None

    return matches, counter
