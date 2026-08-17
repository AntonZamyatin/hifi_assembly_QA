"""Interactive exploration frontend for the assembly-QA report figure.

See `docs/qaview.md` for the design.  This package holds the renderer-agnostic
compute layer (`region`) and, later, the Bokeh/Panel application.  Numbers come from
`qatool` so the interactive view and the published figure can never disagree.
"""

from qaview.region import RegionData

__all__ = ["RegionData"]
