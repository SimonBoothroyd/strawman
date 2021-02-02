import abc
from typing import Any, Dict, Type, TypeVar

import numpy
from pydantic import BaseModel, Field, conint

T = TypeVar("T", bound="WrappedModel")


class IndexRange(BaseModel):
    """Represents a range of indices."""

    start: conint(ge=0) = Field(..., description="The start index.")
    end: conint(ge=1) = Field(..., description="The end index.")

    def __hash__(self):
        return hash((self.start, self.end))

    def __eq__(self, other):
        return (
            type(self) == type(other)
            and self.start == other.start
            and self.end == other.end
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    @classmethod
    def from_key(cls, value: str) -> "IndexRange":
        return IndexRange(start=int(value.split("-")[0]), end=int(value.split("-")[1]))

    def to_key(self) -> str:
        return f"{self.start}-{self.end}"


class ImmutableModel(BaseModel):
    class Config:
        allow_mutation = False
        json_encoders = {numpy.ndarray: lambda x: x.tolist()}


class InnerModel(BaseModel):
    class Config:
        validate_assignment = True
        json_encoders = {numpy.ndarray: lambda x: x.tolist()}


class WrappedModel(abc.ABC):
    class Data(BaseModel):
        pass

    def __init__(self):
        self._inner_data = self.Data()

    def to_json(self) -> str:
        return self._inner_data.json()

    def to_dict(self) -> Dict[str, Any]:
        return self._inner_data.dict()

    @classmethod
    def from_json(cls: Type[T], json_string: str) -> T:
        return_value = cls()
        return_value._inner_data = cls.Data.parse_raw(json_string)

        return return_value

    @classmethod
    def from_dict(cls: Type[T], dict_object: Dict[str, Any]) -> T:
        return_value = cls()
        return_value._inner_data = cls.Data.parse_obj(dict_object)
        return return_value

    def __str__(self):
        return str(self._inner_data)

    def __repr__(self):
        return repr(self._inner_data)
