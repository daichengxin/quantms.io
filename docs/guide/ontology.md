# Ontology Commands

Manage controlled vocabulary (CV) ontology data used by QPX, including PSI-MS and PRIDE CV.

```python exec="1" session="ontology_utils" result="ansi"
import click
import textwrap


def get_click_type_display(param):
    param_type = param.type
    type_str = str(param_type)
    if "Path" in type_str:
        if hasattr(param_type, "dir_okay") and not param_type.dir_okay:
            return "FILE"
        elif hasattr(param_type, "file_okay") and not param_type.file_okay:
            return "DIRECTORY"
        else:
            return "PATH"
    elif isinstance(param_type, click.types.FloatParamType):
        return "FLOAT"
    elif isinstance(param_type, click.types.IntParamType):
        return "INTEGER"
    elif param.is_flag:
        return "FLAG"
    elif isinstance(param_type, click.Choice):
        return "CHOICE"
    else:
        return "TEXT"


def generate_params_table(command):
    table = "<table>\n<thead>\n<tr>\n"
    table += "<th>Parameter</th><th>Type</th><th>Required</th><th>Default</th><th>Description</th>\n"
    table += "</tr>\n</thead>\n<tbody>\n"
    for param in command.params:
        if isinstance(param, click.Option) and param.name not in ["help"]:
            param_names = param.opts
            param_name = param_names[0] if param_names else f"--{param.name}"
            param_type = get_click_type_display(param)
            required = "Yes" if param.required else "No"
            description = param.help or ""
            default_from_help = None
            if "(default:" in description.lower():
                import re

                match = re.search(r"\(default:\s*([^)]+)\)", description, re.IGNORECASE)
                if match:
                    default_from_help = match.group(1).strip()
            if default_from_help:
                default = default_from_help
            elif param.default is not None and str(param.default) != "Sentinel.UNSET":
                if param.is_flag:
                    default = "-"
                elif isinstance(param.default, (int, float)):
                    default = str(param.default)
                elif isinstance(param.default, str):
                    default = f"<code>{param.default}</code>"
                else:
                    default = str(param.default)
            else:
                default = "-"
            table += f"<tr>\n<td><code>{param_name}</code></td>\n<td>{param_type}</td>\n<td>{required}</td>\n<td>{default}</td>\n<td>{description}</td>\n</tr>\n"
    # Also handle Arguments
    for param in command.params:
        if isinstance(param, click.Argument):
            param_name = param.name.upper()
            table += f"<tr>\n<td><code>{param_name}</code></td>\n<td>TEXT</td>\n<td>Yes</td>\n<td>-</td>\n<td>Positional argument</td>\n</tr>\n"
    table += "</tbody>\n</table>"
    return table


def generate_description(command):
    if command.help:
        help_text = command.help
        if "Example" in help_text:
            description = help_text.split("Example")[0].strip()
        else:
            description = help_text.strip()
        lines = description.split("\n")
        if len(lines) > 1:
            description = "\n".join(lines[1:]).strip()
            return f"<p>{description}</p>"
    return ""


def generate_example(command, default_text=""):
    if command.help and "Example" in command.help:
        example_section = command.help.split("Example")[1]
        if ":" in example_section:
            example_section = example_section.split(":", 1)[1]
        example_section = textwrap.dedent(example_section).strip()
        output = ""
        if default_text:
            output += f"<p>{default_text}</p>\n"
        output += f'<pre><code class="language-bash">{example_section}</code></pre>'
        return output
    return ""
```

## Overview

The `ontology` command group manages the controlled vocabulary (CV) ontology data that QPX uses to annotate and validate proteomics metadata. QPX ships with pre-built Parquet files for PSI-MS and PRIDE CV ontologies.

## Available Commands

- [info](#info) - Show loaded ontology sources, versions, and cache status
- [search](#search) - Search for CV terms matching a query string
- [update](#update) - Force download the latest ontology data
- [build](#build) - Rebuild ontology Parquet files from OBO sources (maintainer)

---

## info

Show loaded ontology sources, versions, and cache status.

### Parameters {#info-parameters}

```python exec="1" html="1" session="ontology_utils"
from qpx.cli.ontology import info

print(generate_params_table(info))
```

### Usage Examples {#info-examples}

```python exec="1" html="1" session="ontology_utils"
from qpx.cli.ontology import info

print(generate_example(info, "Show ontology source status:"))
```

---

## search

Search for CV terms matching a query string.

### Parameters {#search-parameters}

```python exec="1" html="1" session="ontology_utils"
from qpx.cli.ontology import search

print(generate_params_table(search))
```

### Usage Examples {#search-examples}

```python exec="1" html="1" session="ontology_utils"
from qpx.cli.ontology import search

print(generate_example(search, "Search for CV terms:"))
```

---

## update

Force download the latest ontology data from the repository.

### Parameters {#update-parameters}

```python exec="1" html="1" session="ontology_utils"
from qpx.cli.ontology import update

print(generate_params_table(update))
```

### Usage Examples {#update-examples}

```python exec="1" html="1" session="ontology_utils"
from qpx.cli.ontology import update

print(generate_example(update, "Update ontology data:"))
```

---

## build

Rebuild ontology Parquet files from OBO sources. This is a maintainer command for regenerating the bundled ontology data.

### Parameters {#build-parameters}

```python exec="1" html="1" session="ontology_utils"
from qpx.cli.ontology import build

print(generate_params_table(build))
```

### Usage Examples {#build-examples}

```python exec="1" html="1" session="ontology_utils"
from qpx.cli.ontology import build

print(generate_example(build, "Rebuild ontology Parquet files:"))
```
