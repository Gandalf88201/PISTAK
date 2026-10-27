"""Registro dei parametri (vedi `registry.py`). Importare questo pacchetto registra quelli predefiniti."""
from pistak.params import builtin as _builtin  # noqa: F401  (registra i parametri del sito)
from pistak.params import extra as _extra  # noqa: F401  (registra i parametri della fase 2)
from pistak.params.plugins import load_errors, load_plugins, plugin_dirs, plugin_files  # noqa: F401
from pistak.params.registry import (  # noqa: F401
    REGISTRY,
    ParamSpec,
    Registry,
    compute_all,
    get_spec,
    parameter,
    sources,
    specs,
)
from pistak.params.structure import (  # noqa: F401
    STRUCTURE_REGISTRY,
    StructureRegistry,
    StructureSpec,
    compute_structure,
    get_structure_spec,
    structure_parameter,
    structure_specs,
)
