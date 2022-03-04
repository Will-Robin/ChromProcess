import numpy as np


class DataReport:
    def __init__(self):
        """
        Parameters
        ----------
        file: str or pathlib Path
        """

        self.filename = "not specified"
        self.experiment_code = "not specified"
        self.conditions = {}
        self.analysis_details = {}
        self.series_values = np.array([])
        self.errors = np.array([])
        self.series_unit = "not specified"
        self.data = {}

    def write_to_file(self, filename=""):
        """
        Write the data report to a csv file.

        Parameters
        ----------
        filename: str
            name for file
        path: pathlib Path object
            Path to folder for file storage.

        Returns
        -------
        None
        """

        import ChromProcess.Writers as write

        write.data_report_to_csv(self, filename=filename)

    def find_repeat_data_entries(self):
        """
        Find data entries that have similar assignments.

        Parameters
        ----------

        Returns
        -------
        repeat_entries: list
            A list of the names of entries which are repeated.
        """

        entries = []
        repeat_entries = []
        for d in self.data:
            token = d.split(" ")[0]
            if token in entries:
                repeat_entries.append(token)
            entries.append(token)

        repeated_entries = list(set(repeat_entries))

        return repeated_entries

    def remove_repeat_entries(self):
        """
        If entries are repeated, remove the one with the lower total magnitude.
        Modified the data report in place.

        TODO: remove? This could be too specific.

        Parameters
        ----------

        Returns
        -------
        None
        """

        # deleting duplicate entries: taking the entry with the higher signal using the
        # signal sum as a discriminant.
        repeat_entries = self.find_repeat_data_entries()

        for r in repeat_entries:
            compare_keys = []
            for d in self.data:
                if r in d:
                    compare_keys.append(d)

            checkline = np.zeros(len(compare_keys))
            for c, comp in enumerate(compare_keys):
                checkline[c] = np.sum(self.data[comp])

            i_min = np.argmin(checkline)

            del self.data[compare_keys[i_min]]

    def remove_specific_entries(self, remove_list):
        """
        Remove the entries in the data report's data using list of them.

        Parameters
        ----------
        remove_list: list
            List of entries to remove from self.data

        Returns
        -------
        None
        """

        for r in remove_list:
            del self.data[r]

    def remove_entries_below_threshold(self, threshold):
        """
        Remove entries of the data report whose maximal values do not exceed a
        threshold.

        Parameters
        ----------
        threshold: float
            threshold below which entries will be removed.

        Returns
        -------
        None
        """

        # remove entries whose concentrations/integrals do not cross a defined boundary
        del_list = []
        for d in self.data:
            if np.amax(self.data[d]) < threshold:
                del_list.append(d)

        self.remove_specific_entries(del_list)
