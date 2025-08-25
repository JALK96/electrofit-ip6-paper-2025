"""Domain logic for computing and applying atomic partial charges.

Refactored functionality now lives in ``process_conformer`` and related helpers. Future enhancements may split
Gaussian/RESP orchestration into a higher-level service object if additional
stateful coordination emerges.

Guidelines:
	* Avoid direct subprocess calls (delegate to adapters layer).
	* Keep I/O side-effects explicit via function parameters.
"""

__all__ = []
