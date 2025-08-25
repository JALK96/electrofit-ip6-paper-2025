"""Curly brace annotation helper (umgezogen aus utils_curly_brace).

Public API: draw_curly_brace (alias für frühere curlyBrace Funktion).
"""
from __future__ import annotations

from typing import Any, Tuple, List
import numpy as np

# Original Funktionen minimal angepasst: PEP8 Name + __all__

__all__ = ["draw_curly_brace"]


def _get_axes_pixel_size(fig, ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    ax_width, ax_height = bbox.width * fig.dpi, bbox.height * fig.dpi
    return ax_width, ax_height


def draw_curly_brace(
    fig,
    ax,
    p1,
    p2,
    k_r: float = 0.1,
    auto: bool = True,
    text: str = "",
    text_offset_lines: int = 2,
    fontdict: dict | None = None,
    **kwargs: Any,
) -> Tuple[float, List[float], list, list, list, list]:
    fontdict = fontdict or {}
    """
    .. _curlyBrace :

    Plot an optionally annotated curly bracket on the given axes of the given figure.

    Note that the brackets are anti-clockwise by default. To reverse the text position, swap
    "p1" and "p2".

    Note that, when the axes aspect is not set to "equal", the axes coordinates need to be
    transformed to screen coordinates, otherwise the arcs may not be seeable.

    Parameters
    ----------
    fig : matplotlib figure object
        The of the target axes.

    ax : matplotlib axes object
        The target axes.

    p1 : two element numeric list
        The coordinates of the starting point.

    p2 : two element numeric list
        The coordinates of the end point.

    k_r : float
        This is the gain controlling how "curvy" and "pointy" (height) the bracket is.

        Note that, if this gain is too big, the bracket would be very strange.

    bool_auto : boolean
        This is a switch controlling wether to use the auto calculation of axes
        scales.

        When the two axes do not have the same aspects, i.e., not "equal" scales,
        this should be turned on, i.e., True.

        When "equal" aspect is used, this should be turned off, i.e., False.

        If you do not set this to False when setting the axes aspect to "equal",
        the bracket will be in funny shape.

        Default = True

    str_text : string
        The annotation text of the bracket. It would displayed at the mid point
        of bracket with the same rotation as the bracket.

        By default, it follows the anti-clockwise convention. To flip it, swap
        the end point and the starting point.

        The appearance of this string can be set by using "fontdict", which follows
        the same syntax as the normal matplotlib syntax for font dictionary.

        Default = empty string (no annotation)

    int_line_num : int
        This argument determines how many lines the string annotation is from the summit
        of the bracket.

        The distance would be affected by the font size, since it basically just a number of
        lines appended to the given string.

        Default = 2

    fontdict : dictionary
        This is font dictionary setting the string annotation. It is the same as normal
        matplotlib font dictionary.

        Default = empty dict

    **kwargs : matplotlib line setting arguments
        This allows the user to set the line arguments using named arguments that are
        the same as in matplotlib.

    Returns
    -------
    theta : float
        The bracket angle in radians.

    summit : list
        The positions of the bracket summit.

    arc1 : list of lists
        arc1 positions.

    arc2 : list of lists
        arc2 positions.

    arc3 : list of lists
        arc3 positions.

    arc4 : list of lists
        arc4 positions.

    Reference
    ----------
    https://uk.mathworks.com/matlabcentral/fileexchange/38716-curly-brace-annotation
    """
    pt1 = [None, None]
    pt2 = [None, None]

    ax_width, ax_height = _get_axes_pixel_size(fig, ax)
    ax_xlim = list(ax.get_xlim())
    ax_ylim = list(ax.get_ylim())

    # log scale x
    if "log" in ax.get_xaxis().get_scale():
        pt1[0] = np.log(p1[0]) if p1[0] > 0 else (-np.log(abs(p1[0])) if p1[0] < 0 else 0.0)
        pt2[0] = np.log(p2[0]) if p2[0] > 0 else (-np.log(abs(p2[0])) if p2[0] < 0 else 0.0)
        for i, v in enumerate(ax_xlim):
            ax_xlim[i] = np.log(v) if v > 0 else (-np.log(abs(v)) if v < 0 else 0.0)
    else:
        pt1[0], pt2[0] = p1[0], p2[0]

    # log scale y
    if "log" in ax.get_yaxis().get_scale():
        pt1[1] = np.log(p1[1]) if p1[1] > 0 else (-np.log(abs(p1[1])) if p1[1] < 0 else 0.0)
        pt2[1] = np.log(p2[1]) if p2[1] > 0 else (-np.log(abs(p2[1])) if p2[1] < 0 else 0.0)
        for i, v in enumerate(ax_ylim):
            ax_ylim[i] = np.log(v) if v > 0 else (-np.log(abs(v)) if v < 0 else 0.0)
    else:
        pt1[1], pt2[1] = p1[1], p2[1]

    xscale = ax_width / abs(ax_xlim[1] - ax_xlim[0])
    yscale = ax_height / abs(ax_ylim[1] - ax_ylim[0])
    if not auto:
        xscale = yscale = 1.0

    pt1[0] = (pt1[0] - ax_xlim[0]) * xscale
    pt1[1] = (pt1[1] - ax_ylim[0]) * yscale
    pt2[0] = (pt2[0] - ax_xlim[0]) * xscale
    pt2[1] = (pt2[1] - ax_ylim[0]) * yscale

    theta = np.arctan2(pt2[1] - pt1[1], pt2[0] - pt1[0])
    r = np.hypot(pt2[0] - pt1[0], pt2[1] - pt1[1]) * k_r

    x11 = pt1[0] + r * np.cos(theta)
    y11 = pt1[1] + r * np.sin(theta)
    x22 = (pt2[0] + pt1[0]) / 2.0 - 2.0 * r * np.sin(theta) - r * np.cos(theta)
    y22 = (pt2[1] + pt1[1]) / 2.0 + 2.0 * r * np.cos(theta) - r * np.sin(theta)
    x33 = (pt2[0] + pt1[0]) / 2.0 - 2.0 * r * np.sin(theta) + r * np.cos(theta)
    y33 = (pt2[1] + pt1[1]) / 2.0 + 2.0 * r * np.cos(theta) + r * np.sin(theta)
    x44 = pt2[0] - r * np.cos(theta)
    y44 = pt2[1] - r * np.sin(theta)

    q = np.linspace(theta, theta + np.pi / 2.0, 50)
    t = q[::-1]

    arc1x = r * np.cos(t + np.pi / 2.0) + x11
    arc1y = r * np.sin(t + np.pi / 2.0) + y11
    arc2x = r * np.cos(q - np.pi / 2.0) + x22
    arc2y = r * np.sin(q - np.pi / 2.0) + y22
    arc3x = r * np.cos(q + np.pi) + x33
    arc3y = r * np.sin(q + np.pi) + y33
    arc4x = r * np.cos(t) + x44
    arc4y = r * np.sin(t) + y44

    arc1x = arc1x / xscale + ax_xlim[0]
    arc2x = arc2x / xscale + ax_xlim[0]
    arc3x = arc3x / xscale + ax_xlim[0]
    arc4x = arc4x / xscale + ax_xlim[0]
    arc1y = arc1y / yscale + ax_ylim[0]
    arc2y = arc2y / yscale + ax_ylim[0]
    arc3y = arc3y / yscale + ax_ylim[0]
    arc4y = arc4y / yscale + ax_ylim[0]

    if "log" in ax.get_xaxis().get_scale():
        for arr in (arc1x, arc2x, arc3x, arc4x):
            for i, vx in enumerate(arr):
                arr[i] = np.exp(vx) if vx > 0 else (-np.exp(abs(vx)) if vx < 0 else 0.0)
    if "log" in ax.get_yaxis().get_scale():
        for arr in (arc1y, arc2y, arc3y, arc4y):
            for i, vy in enumerate(arr):
                arr[i] = np.exp(vy) if vy > 0 else (-np.exp(abs(vy)) if vy < 0 else 0.0)

    ax.plot(arc1x, arc1y, **kwargs)
    ax.plot(arc2x, arc2y, **kwargs)
    ax.plot(arc3x, arc3y, **kwargs)
    ax.plot(arc4x, arc4y, **kwargs)
    ax.plot([arc1x[-1], arc2x[1]], [arc1y[-1], arc2y[1]], **kwargs)
    ax.plot([arc3x[-1], arc4x[1]], [arc3y[-1], arc4y[1]], **kwargs)

    summit = [arc2x[-1], arc2y[-1]]

    if text:
        offset_lines = int(text_offset_lines)
        spacer = "\n" * offset_lines
        ang_deg = float(np.degrees(theta) % 360.0)
        if 0.0 <= ang_deg <= 90.0:
            rotation = ang_deg
            text_out = text + spacer
        elif 90.0 < ang_deg < 270.0:
            rotation = ang_deg + 180.0
            text_out = spacer + text
        else:  # 270..360
            rotation = ang_deg
            text_out = text + spacer
        ax.axes.text(
            arc2x[-1],
            arc2y[-1],
            text_out,
            ha="center",
            va="center",
            rotation=rotation,
            fontdict=fontdict,
        )

    return theta, summit, [arc1x, arc1y], [arc2x, arc2y], [arc3x, arc3y], [arc4x, arc4y]