from typing import Literal

from pydantic import conint

AtomicNumber = conint(ge=1, le=118)
AtomicSymbol = Literal["H", "C", "O", "N"]


def atomic_number_to_symbol(atomic_number: int) -> AtomicSymbol:
    return {1: "H", 6: "C", 7: "N", 8: "O"}[atomic_number]
