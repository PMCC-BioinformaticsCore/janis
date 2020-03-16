from typing import List, Dict
from datetime import date
from distutils.version import StrictVersion

from janis_core import Tool, Metadata
from requests.utils import requote_uri


class TocObject:
    def __init__(self, title, description, url):
        self.title = title
        self.description = description
        self.url = url


class NestedDictionaryTypeException(Exception):
    def __init__(self, key, error_type, keys=None):
        keychain = (" for keys: " + ".".join(keys)) if keys else ""
        super(NestedDictionaryTypeException, self).__init__(
            f"Incorrect type '{error_type}' for key '{key}'{keychain}"
        )
        self.keys = keys
        self.key = key
        self.error_type = error_type


def format_rst_link(text, link):
    return f"`{text} <{link}>`_"


def get_tool_url(toolname, version):
    if not (toolname and version):
        return ""
    return requote_uri(toolname + "_" + str(version)).lower()


def prepare_byline(
    shortdocumentation: str, contributors: List[str], versions: List[str]
):
    short_prepared = (
        shortdocumentation[:-1]
        if shortdocumentation and shortdocumentation[-1] == "."
        else shortdocumentation
    )

    contributors, versions = contributors or [], versions or []
    contributorstr = (
        f"{len(contributors)} contributor{'s' if len(contributors) != 1 else ''}"
    )
    versionstr = f"{len(versions)} version{'s' if len(versions) != 1 else ''}"
    return (
        "*"
        + " · ".join(t for t in [short_prepared, contributorstr, versionstr] if t)
        + "*"
    )


def nested_keys_append_with_root(d, keys: List[str], value, root_key):
    if len(keys) == 0:
        if root_key in d:
            d[root_key].append(value)
        else:
            d[root_key] = [value]
        return d
    try:
        key = keys[0]
        if key in d:
            if not isinstance(d[key], dict):
                raise NestedDictionaryTypeException(
                    key=key, error_type=type(d[key]), keys=keys
                )
            nested_keys_append_with_root(d[key], keys[1:], value, root_key=root_key)
            return d
        else:
            d[key] = nested_keys_append_with_root(
                {}, keys[1:], value, root_key=root_key
            )
    except NestedDictionaryTypeException as de:
        raise NestedDictionaryTypeException(
            key=de.key, error_type=de.error_type, keys=keys
        )

    return d


def nested_keys_add(d, keys: List[str], value) -> bool:
    if len(keys) == 0:
        raise Exception(
            f"Couldn't add {value} to nested dictionary as no keys were provided"
        )
    key = keys[0]
    if len(keys) == 1:
        if key in d:
            return False
        d[key] = value
        return True
    try:
        if key not in d:
            d[key] = {}
        elif not isinstance(d[key], dict):
            raise NestedDictionaryTypeException(
                key=key, error_type=type(d[key]), keys=keys
            )
        return nested_keys_add(d[key], keys[1:], value)
    except NestedDictionaryTypeException as de:
        raise NestedDictionaryTypeException(
            key=de.key, error_type=de.error_type, keys=keys
        )


def get_toc(title, intro_text, subpages, caption="Contents", max_depth=1):
    prepared_subpages = "\n".join(
        "   " + m.lower() for m in sorted(subpages, key=lambda l: l.lower())
    )
    return f"""
{title.replace('{title}', title)}
{"=" * len(title)}

{intro_text}

.. toctree::
   :maxdepth: {max_depth}
   :caption: {caption}:

{prepared_subpages}

*This page was auto-generated on {date.today().strftime(
        "%d/%m/%Y")}. Please do not directly alter the contents of this page.*
"""


def get_tool_toc(
    alltoolsmap: Dict[str, Dict[str, Tool]],
    title,
    intro_text,
    subpages,
    tools,
    caption="Contents",
    max_depth=1,
):

    pd = " " * 5
    mappedtools = "\n".join(
        "\n".join(pd + r for r in get_tool_row(alltoolsmap[tool]).split("\n"))
        for tool in tools
    )

    mappedmodules = "\n\n".join(
        f'     <li><a href="{m.lower()}/index.html">{m}</a></li>'
        for m in sorted(subpages)
    )

    return f"""
:orphan:

{title.replace('{title}', title)}
{"=" * len(title)}

{intro_text}

.. raw:: html

   <ul>
{mappedmodules}
   </ul>
{mappedtools}

"""


def get_tool_row(tools: Dict[str, Tool]):
    distincted = {t.version(): t for t in tools.values()}
    versions = sort_tool_versions(list(distincted.keys()))
    latestversion = versions[0]

    tool = distincted[latestversion]
    meta: Metadata = tool.bind_metadata() or tool.metadata
    sd = meta.short_documentation
    sdstr = f'<p style="color: black; margin-bottom: 10px">{sd}' if sd else ""

    href = tool.id().lower() + ".html"
    return f"""\
<a href="{href}">
  <p style="margin-bottom: 5px"><b>{tool.friendly_name()}</b> <span style="margin-left: 10px; color: darkgray">{tool.id()}</span></p>
  {sdstr}
  <p><span style="margin-right: 10px; color: darkgray">({len(versions)} versions)</span>{version_html(latestversion, href=href)}</p>
</a>
<hr />
    """


def sort_tool_versions(versions: List[str]) -> List[str]:
    try:
        return sorted(versions, key=StrictVersion, reverse=True)
    except:
        return sorted(versions, reverse=True)


def version_html(version, href=None):
    fhref = f'href="{href}"' if href else ""
    return f"""\
<a class="version-button" {fhref} style="margin-bottom: 10px">
  v<b>{version[1:] if (version and version.startswith("v")) else version}</b>
</a>"""
