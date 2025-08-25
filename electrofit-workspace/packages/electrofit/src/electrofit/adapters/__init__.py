"""Adapters für externe Tools (GROMACS etc.)."""
from .gromacs import set_up_production  # re-export primäre Funktion
__all__ = ["set_up_production"]
