from typing import Literal, Any
from pydantic import BaseModel, Field, field_serializer
from .records import BaseRecord
from .analysis import BaseAnalysis

class Manifest(BaseModel):
    hich_version: Literal["unstable"] = "unstable"
    records: dict[str, BaseRecord] = {}
    analysis: dict[str, BaseAnalysis] = {}
    __template: BaseRecord = None

    def set_template(self, template: BaseRecord):
        self.__template = template
    
    def add_record(self, record: dict | BaseRecord):
        if isinstance(record, BaseRecord):
            id = record.id
        elif self.__template is not None:
            record_dumped = (
                self.__template
                .model_copy(update=record)
                .model_dump()
            )
            RECORD_TYPE = type(self.__template)
            record = RECORD_TYPE.model_validate(record_dumped)
            id = record.id
        self.records[id] = record
    
    def add_analysis(self, analysis: BaseAnalysis):
        if isinstance(analysis, BaseAnalysis):
            analysis_type = analysis.analysis_type
            self.analysis.setdefault(analysis_type, [])
            self.analysis[analysis_type].append(analysis)

    @field_serializer('records')
    def serialize_records(self, records: dict[str, BaseRecord]) -> dict[str, Any]:
        """
        Forces Pydantic to call model_dump() on each *instance*,
        thus using its specific subclass fields.
        """
        # We must return a dict, but the values can be 'Any'
        # (which will be the result of each record's own model_dump)
        return {key: record.model_dump() for key, record in records.items()}
    
    @field_serializer('analysis')
    def serialize_analysis(self, all_analysis: dict[str, BaseAnalysis]) -> dict[str, Any]:
        """
        Forces Pydantic to call model_dump() on each *instance*,
        thus using its specific subclass fields.
        """
        # We must return a dict, but the values can be 'Any'
        # (which will be the result of each record's own model_dump)
        serialized = {}
        for analysis_type, analysis_for_type in all_analysis.items():
            serialized[analysis_type] = []
            for analysis in analysis_for_type:
                serialized[analysis_type].append(analysis.model_dump())

        return serialized