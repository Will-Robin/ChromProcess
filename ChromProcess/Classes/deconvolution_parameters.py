class Deconvolution:
    def __init__(
        self, parameters: dict[str, float], region: list[float], function: str
    ):
        self.parameters: dict[str, float] = parameters
        self.region: list[float] = region
        self.function: str = function
