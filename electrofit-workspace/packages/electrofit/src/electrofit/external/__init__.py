"""Deprecated namespace: will be renamed to electrofit.adapters.

Existing imports continue to work during deprecation window.
"""
import warnings as _warnings
_warnings.warn(
	"Namespace 'electrofit.external' wird zu 'electrofit.adapters' umbenannt.",
	DeprecationWarning,
	stacklevel=2,
)
