"""Plot helpers for QPX views — lazy plotly import."""

from __future__ import annotations


def _get_plotly():
    """Lazy plotly import."""
    try:
        import plotly.graph_objects as go

        return go
    except ImportError:
        raise ImportError("plotly is required for plotting. Install it with: pip install qpx[plotting]")


def bar_chart(
    df,
    x: str,
    y: str,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    rotation: int = 45,
):
    """Simple bar chart from a DataFrame."""
    go = _get_plotly()
    fig = go.Figure(go.Bar(x=df[x].astype(str), y=df[y]))
    fig.update_layout(
        title=title,
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        xaxis_tickangle=-rotation,
    )
    return fig


def scatter_plot(
    df,
    x: str,
    y: str,
    color_by: str | None = None,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
):
    """Scatter plot from a DataFrame, optionally colored by a grouping column."""
    go = _get_plotly()
    fig = go.Figure()
    if color_by and color_by in df.columns:
        for group_val in df[color_by].unique():
            mask = df[color_by] == group_val
            fig.add_trace(
                go.Scatter(
                    x=df.loc[mask, x],
                    y=df.loc[mask, y],
                    mode="markers",
                    name=str(group_val),
                )
            )
    else:
        fig.add_trace(go.Scatter(x=df[x], y=df[y], mode="markers"))
    fig.update_layout(title=title, xaxis_title=xlabel, yaxis_title=ylabel)
    return fig


def grouped_bar_chart(
    df,
    x: str,
    y_cols: list[str],
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
):
    """Grouped bar chart with multiple y-columns."""
    go = _get_plotly()
    fig = go.Figure()
    x_vals = df[x].astype(str)
    for col in y_cols:
        fig.add_trace(go.Bar(x=x_vals, y=df[col], name=col))
    fig.update_layout(
        barmode="group",
        title=title,
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        xaxis_tickangle=-45,
    )
    return fig
