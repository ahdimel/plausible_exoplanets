"""Plausibility flagging framework.

Every generated object carries a list of Flag records. Severities:

- INFO         : notable property, fully consistent with known physics/data
- QUESTIONABLE : allowed by physics but in a regime where models or data are
                 weak; the object is kept and the doubt is recorded
- INVALID      : violates physics or robust empirical limits; the generator
                 must reject and resample (INVALID flags never appear in the
                 output database except in generator diagnostics)
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import List


class Severity(str, Enum):
    INFO = "info"
    QUESTIONABLE = "questionable"
    INVALID = "invalid"


@dataclass
class Flag:
    severity: Severity
    rule: str          # machine-readable rule id, e.g. "density.exceeds_pure_iron"
    message: str       # human-readable explanation

    def to_dict(self) -> dict:
        d = asdict(self)
        d["severity"] = self.severity.value
        return d


@dataclass
class Flagged:
    """Mixin-style container for objects that accumulate flags."""

    flags: List[Flag] = field(default_factory=list)

    def add_flag(self, severity: Severity, rule: str, message: str) -> None:
        self.flags.append(Flag(severity, rule, message))

    @property
    def is_invalid(self) -> bool:
        return any(f.severity == Severity.INVALID for f in self.flags)

    @property
    def questionable_flags(self) -> List[Flag]:
        return [f for f in self.flags if f.severity == Severity.QUESTIONABLE]
