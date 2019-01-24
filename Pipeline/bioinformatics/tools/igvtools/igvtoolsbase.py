from abc import abstractmethod, ABC


class IgvToolsBase(ABC):

    @staticmethod
    @abstractmethod
    def docker():
        raise Exception("An error likely occurred when resolving the method order for docker for the igvtools classes "
                        "or you're trying to execute the docker method of the base class (ie, don't do that). "
                        "The method order resolution must preference BwaBase subclasses, "
                        "and the subclass must contain a definition for docker.")