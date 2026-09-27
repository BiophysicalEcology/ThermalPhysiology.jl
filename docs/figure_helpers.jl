module FigureHelpers

using CairoMakie, Markdown, Unitful

export figure_axis, parameter_table, markdown_table

"""
    figure_axis(xlabel, ylabel; size=(700, 500), kw...)

A `Figure` and `Axis` with minor ticks and grid lines, as used for the figures of the manual.
"""
function figure_axis(xlabel, ylabel; size=(700, 500), kw...)
    fig = Figure(; size)
    ax = Axis(fig[1, 1];
        xlabel, ylabel,
        xminorticksvisible=true, yminorticksvisible=true,
        xminorgridvisible=true, yminorgridvisible=true,
        xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
        kw...,
    )
    return fig, ax
end

"""
    parameter_table(model)
    parameter_table(label => model, ...)

A Markdown table of the parameters (fields) of one model, or of several models of the same type side by side,
one column per `label => model`, with values rounded to 4 significant digits.
"""
parameter_table(model) = parameter_table("Value" => model)
function parameter_table(columns::Pair...)
    models = last.(columns)
    header = "| Parameter | " * join(first.(columns), " | ") * " |"
    separator = "|:--|" * repeat(":--|", length(columns))
    rows = ["| `$name` | " * join((_format(getfield(model, name)) for model in models), " | ") * " |"
            for name in fieldnames(typeof(first(models)))]
    return Markdown.parse(join([header, separator, rows...], "\n"))
end

"""
    markdown_table(header, rows)

A Markdown table with column names `header` and one row per element of `rows`, with numbers and quantities
rounded to 4 significant digits.
"""
function markdown_table(header, rows)
    lines = ["| " * join(header, " | ") * " |", "|" * repeat(":--|", length(header))]
    append!(lines, ["| " * join(map(_format, row), " | ") * " |" for row in rows])
    return Markdown.parse(join(lines, "\n"))
end

_format(x::Unitful.Quantity) = string(round(unit(x), x; sigdigits = 4))
_format(x::Real) = string(round(x; sigdigits = 4))
_format(::Nothing) = "–"
_format(x::AbstractMatrix) = "$(join(size(x), " × ")) matrix"
_format(x) = string(x)

end
