from . import cutout
from . import depth

def __getattr__(name):
    """Get realizations using lazy import from
    `PEP 562 <https://www.python.org/dev/peps/pep-0562/>`_.
    Raises
    ------
    AttributeError
        If "name" is not in :mod:`astropy.cosmology.realizations`
    """
    if name not in realizations.available:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}.")

    return getattr(realizations, name)


def __dir__():
    """Directory, including lazily-imported objects."""
    return 
