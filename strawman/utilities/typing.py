from typing import TYPE_CHECKING, Any, get_args

import numpy
from openff.units import unit


class _FloatQuantityMeta(type):
    def __getitem__(self, t):
        return type("FloatQuantity", (FloatQuantity,), {"__unit__": t})


class FloatQuantity(float, metaclass=_FloatQuantityMeta):
    @classmethod
    def __get_validators__(cls):
        yield cls.validate_type

    @classmethod
    def validate_type(cls, val):

        if isinstance(val, float):
            return val

        expected_unit = unit.Unit(getattr(cls, "__unit__", Any))

        if isinstance(val, unit.Quantity):
            return val.to(expected_unit)

        raise NotImplementedError()


class ArrayQuantityMeta(type):
    def __getitem__(self, t):
        return type("ArrayQuantity", (ArrayQuantity,), {"__unit__": t})


class ArrayQuantity(float, metaclass=ArrayQuantityMeta):
    @classmethod
    def __get_validators__(cls):
        yield cls.validate_type

    @classmethod
    def validate_type(cls, val):

        if isinstance(val, list):
            val = numpy.array(val)

        if isinstance(val, numpy.ndarray):
            return val

        unit_type = get_args(getattr(cls, "__unit__"))[0]
        expected_unit = unit.Unit(unit_type)

        if isinstance(val, unit.Quantity):
            val = val.to(expected_unit).magnitude
            assert isinstance(
                val, numpy.ndarray
            ), f"invalid inner type of {type(val)}, expected a numpy array."

            return val

        raise NotImplementedError()


if TYPE_CHECKING:
    FloatQuantity = float
    ArrayQuantity = numpy.ndarray
