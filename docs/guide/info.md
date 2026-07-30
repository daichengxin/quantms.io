# Info Commands

Inspect QPX datasets — view summaries, schemas, and Parquet metadata.

```python exec="1" session="info_utils" result="ansi"
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

The `info` command group provides tools for inspecting QPX datasets. When invoked without a subcommand, it displays a summary of the dataset including available structures and row counts. Subcommands provide detailed schema and Parquet metadata inspection.

## Available Commands

- [info](#info) (default) - Show dataset summary with structures and row counts
- [schema](#schema) - Show Arrow schema for a data structure
- [metadata](#metadata) - Show Parquet footer metadata

---

## info

Show a summary of a QPX dataset.

### Parameters {#info-parameters}

```python exec="1" html="1" session="info_utils"
from qpx.cli.info import info

print(generate_params_table(info))
```

### Description {#info-description}

```python exec="1" html="1" session="info_utils"
from qpx.cli.info import info

print(generate_description(info))
```

### Usage Examples {#info-examples}

```python exec="1" html="1" session="info_utils"
from qpx.cli.info import info

print(generate_example(info, "Show dataset summary:"))
```

---

## schema

Show the Arrow schema for a QPX data structure.

### Description {#schema-description}

```python exec="1" html="1" session="info_utils"
from qpx.cli.info import info_schema_cmd

print(generate_description(info_schema_cmd))
```

### Parameters {#schema-parameters}

```python exec="1" html="1" session="info_utils"
from qpx.cli.info import info_schema_cmd

print(generate_params_table(info_schema_cmd))
```

### Usage Examples {#schema-examples}

```python exec="1" html="1" session="info_utils"
from qpx.cli.info import info_schema_cmd

print(generate_example(info_schema_cmd, "Inspect data structure schemas:"))
```

---

## metadata

Show Parquet footer metadata for QPX files.

### Description {#metadata-description}

```python exec="1" html="1" session="info_utils"
from qpx.cli.info import info_metadata_cmd

print(generate_description(info_metadata_cmd))
```

### Parameters {#metadata-parameters}

```python exec="1" html="1" session="info_utils"
from qpx.cli.info import info_metadata_cmd

print(generate_params_table(info_metadata_cmd))
```

### Usage Examples {#metadata-examples}

```python exec="1" html="1" session="info_utils"
from qpx.cli.info import info_metadata_cmd

print(generate_example(info_metadata_cmd, "View Parquet file metadata:"))
```

---

## Best Practices

- Use `qpxc info` to get a quick overview of a dataset before querying
- Use `qpxc info schema --canonical` to compare the on-disk schema against the QPX specification
- Use `qpxc info metadata` to check Parquet compression, row groups, and key-value metadata
