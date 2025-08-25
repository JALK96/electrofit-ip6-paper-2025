"""Deprecated shim: use electrofit.adapters.gromacs statt electrofit.external.gromacs.

Beibehaltung für rückwärtskompatible CLI-Einstiegspunkte.
"""
from electrofit.adapters.gromacs import set_up_production  # re-export

__all__ = ["set_up_production"]
