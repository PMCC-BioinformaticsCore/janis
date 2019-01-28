import json
from _sha1 import sha1
from datetime import date


class Metadata(object):

    DATE_FORMAT = "%Y-%m-%d"

    def __init__(self, creator=None, maintainer=None, maintainer_email=None, date_created=None, date_updated=None,
                 institution=None, doi=None, citation=None, keywords=None, documentation_url=None, documentation=None,
                 version=None):
        """

        :param creator:
        :param maintainer:
        :param maintainer_email:
        :param date_created:
        :param date_updated:
        :param institution:
        :param doi:
        :param citation:
        :type citation: str | list[str] | None
        :param keywords:
        :param documentation_url:
        :param documentation:
        :param version:
        """
        self.creator = creator
        self.maintainer = maintainer
        self.maintainerEmail = maintainer_email
        self.dateCreated = date_created.strftime(self.DATE_FORMAT) if isinstance(date_created, date) else date_created
        self.dateUpdated = date_created.strftime(self.DATE_FORMAT) if isinstance(date_updated, date) else date_updated
        self.institution = institution
        self.doi = doi
        self.citation = citation
        self.keywords = keywords
        self.documentation = documentation
        self.documentationUrl = documentation_url
        self.version = version

    def update(self, **kwargs):
        for k in kwargs:
            self.__setattr__(k, kwargs[k])
        return self

    def get_dict(self, object_to_checksum):

        checksum = sha1(json.dumps(object_to_checksum, sort_keys=True).encode('utf-8')).hexdigest()
        # https://stackoverflow.com/q/5884066

        d = {k: v for k, v in vars(self).items() if v is not None}
        d["checksum"] = checksum
        d["dateGenerated"] = date.today().strftime(self.DATE_FORMAT)

        return d


class WorkflowMetadata(Metadata):
    def get_dict(self, object_to_checksum, ninputs=None, nsteps=None, noutputs=None):
        d = super(WorkflowMetadata, self).get_dict(object_to_checksum)

        if ninputs:
            d["numberOfInputs"] = ninputs
        if nsteps:
            d["numberOfSteps"] = nsteps
        if noutputs:
            d["numberOfOutputs"] = noutputs

        return d


class ToolMetadata(Metadata):
    pass

