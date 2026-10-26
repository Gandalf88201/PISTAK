"""Registro dei parametri (vedi `registry.py`). Importare questo pacchetto registra quelli predefiniti."""
from pistak.params import builtin as _builtin  # noqa: F401  (registra i parametri del sito)
from pistak.params import extra as _extra  # noqa: F401  (registra i parametri della fase 2)
from pistak.params.plugins import load_errors, load_plugins, plugin_files  # noqa: F401
from pistak.params.registry import ParamSpec, compute_all, errors, get_spec, parameter, specs  # noqa: F401
