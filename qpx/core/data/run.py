"""Run data structure — MS acquisition run metadata."""

from qpx.core.data.base import BaseStructure
from qpx.core.data.loader import load_schema

RunSchema = load_schema("run")


class Run(BaseStructure):
    """MS acquisition run metadata."""

    _schema_class = RunSchema
