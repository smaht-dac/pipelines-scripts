import json
import zipfile
from typing import List

########################################################################
# QMGeneric
# QMValue
########################################################################


class QMValue:
    """
    """

    def __init__(
        self,
        key,
        value,
        visible=None,
        derived_from=None,
        flag=None,
        tooltip=None
    ):
        # Required
        self.key = key
        self.value = value
        # Optional
        self.visible = visible
        self.derived_from = derived_from
        self.flag = flag
        self.tooltip = tooltip

    def to_dict(self):
        return {k: v for (k, v) in self.__dict__.items() if v != None}


class QMGeneric:
    """
    """

    def __init__(self, name):
        self.name = name
        # List of QMValue objects
        self.qm_values = []

    def add_value(self, QMValue_: QMValue):
        self.qm_values.append(QMValue_)

    def add_values(self, QMValues: List[QMValue]):
        self.qm_values += QMValues

    def to_dict(self):
        return {
            'name': self.name,
            'qc_values': [QMValue_.to_dict() for QMValue_ in self.qm_values]
        }

    def create_archive(self, filenames, archive_name):
        with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            for filename in filenames:
                archive.write(filename)

    def write_json(self, filename, dict=None):
        if not dict:
            dict = self.to_dict()
        with open(filename, 'w') as fo:
            json.dump(dict, fo, sort_keys=True, indent=2)
