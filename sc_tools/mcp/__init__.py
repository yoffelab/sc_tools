"""MCP (Model Context Protocol) servers for sc_tools.

Two servers are provided:

``sc-tools``
    Analysis tools: validate checkpoints, generate QC reports,
    score gene signatures, generate SLURM scripts, etc.

``sc-registry``
    Registry bookkeeping: query datasets, SLURM jobs, agent tasks.

Start servers
-------------
::

    # sc-tools analysis tools
    python -m sc_tools.mcp.tools_server

    # sc-registry bookkeeping
    python -m sc_tools.mcp.registry_server

Configure in ``.mcp.json`` at repo root so Claude Code discovers them
automatically.
"""
