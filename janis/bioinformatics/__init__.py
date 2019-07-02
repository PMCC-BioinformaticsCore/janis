try:
    from janis_bioinformatics import tools, data_types
except ModuleNotFoundError:
    print(
        "You must install the bioinformatics tools to use this section:\n\t`janis-pipelines[bioinformatics]` OR\n\t`janis-pipelines.bioinformatics`"
    )
