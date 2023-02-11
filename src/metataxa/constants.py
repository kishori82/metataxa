from enum import Enum
from typing import Type


class DeltaType(Enum):
    DELTA = "DELTA"
    DELTA_STAR = "DELTA_STAR"
    DELTA_PLUS = "DELTA_PLUS"
