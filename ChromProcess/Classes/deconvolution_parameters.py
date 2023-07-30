class Deconvolution:
    def __init__(
        self, parameters: dict[str, float], region: list[float], function: str
    ):
        self.parameters: dict[str, float] = parameters
        self.region: list[float] = region
        self.function: str = function

    def __str__(self):
        parameters = [f"{k}: {self.parameters[k]}" for k in self.parameters]
        str_repr = f"""
        function: {self.function}
        region: {self.region[0]} - {self.region[1]}
        """
        str_repr += "\n".join(parameters)
        return str_repr
