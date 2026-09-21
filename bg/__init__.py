from importlib import import_module
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
	from .BG2CellML import BG2CellMLV1, write_cellmlV1
	from .SCAP import SCAPEngine
	from .defineBG import (
		BondGraph,
		ComponentType,
		ConnectionType,
		Domain,
		DomainRefiner,
		exportBG,
		importBG,
	)
	from .mathBG import EquationBuilder
	from .visualiseBG import drawBG, print_bond_table

__all__ = [
	"BondGraph",
	"ComponentType",
	"Domain",
	"ConnectionType",
	"DomainRefiner",
	"exportBG",
	"importBG",
	"drawBG",
	"print_bond_table",
	"SCAPEngine",
	"EquationBuilder",
	"BG2CellMLV1",
	"write_cellmlV1",
]

_LAZY_EXPORTS = {
	"BondGraph": (".defineBG", "BondGraph"),
	"ComponentType": (".defineBG", "ComponentType"),
    "Domain": (".defineBG", "Domain"),
    "ConnectionType": (".defineBG", "ConnectionType"),
    "DomainRefiner": (".defineBG", "DomainRefiner"),
    "exportBG": (".defineBG", "exportBG"),
    "importBG": (".defineBG", "importBG"),
	"drawBG": (".visualiseBG", "drawBG"),
	"print_bond_table": (".visualiseBG", "print_bond_table"),
	"SCAPEngine": (".SCAP", "SCAPEngine"),
	"EquationBuilder": (".mathBG", "EquationBuilder"),
	"BG2CellMLV1": (".BG2CellML", "BG2CellMLV1"),
	"write_cellmlV1": (".BG2CellML", "write_cellmlV1"),
}


def __getattr__(name: str) -> Any:
	"""Load optional public modules only when one of their exports is requested."""
	try:
		module_name, attribute_name = _LAZY_EXPORTS[name]
	except KeyError as error:
		raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from error
	return getattr(import_module(module_name, __name__), attribute_name)