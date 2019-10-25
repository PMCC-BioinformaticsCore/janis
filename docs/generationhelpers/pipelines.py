from datetime import datetime

from janis_core import WorkflowMetadata, Workflow

from docs.generationhelpers.utils import version_html


def generate_pipeline_box(workflow: Workflow, leading_space=""):

    meta: WorkflowMetadata = workflow.bind_metadata() or workflow.metadata

    tag_component = lambda tag: f'<span class="no-select tagelement">{tag}</span>'
    href = f"{workflow.id().lower()}.html"

    tags = "".join(tag_component(t) for t in (meta.keywords or []))
    date: datetime = meta.dateUpdated or meta.dateCreated

    contributors = meta.contributors or []
    max_display_conts = 5
    if len(contributors) < max_display_conts:
        contributorstr = ", ".join(meta.contributors or ["None"])
    else:
        nothers = len(contributors) - max_display_conts + 2
        contributorstr = (
            ", ".join(meta.contributors[: max_display_conts - 2])
            + f" and {nothers} others"
        )

    return "\n".join(
        leading_space + l
        for l in f"""
<div class="col-6" style="margin: 10px; padding: 20px; border: 1px solid #e3e3e3; border-radius: 5px;">
    <h4 style="margin-bottom: 10px"><a href="{href}">{workflow.friendly_name()}</a></h4>
    {f'<p>{tags}</p>' if tags else ""}
    <p>{meta.short_documentation or "<em>Short documentation required</em>"}</p>
    {version_html(workflow.version() or meta.version)}
    <p style="margin-bottom: 0px; font-size: 12px">Contributors: {contributorstr}</p>
</div>""".split(
            "\n"
        )
    )
