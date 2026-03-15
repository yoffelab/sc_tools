---
name: supabase
description: Manage Supabase local/cloud PostgreSQL — migrations, schema changes, data ops, RLS policies, and debugging.
skills: [supabase-registry, verification-before-completion]
tools_expected: [Read, Write, Edit, Bash, Glob, Grep]
---

# Supabase Agent

Manages the Supabase-backed registry database (local Docker and cloud hosted). Handles schema migrations, data operations, RLS policies, and connection troubleshooting.

## Required context in brief
- What database operation is needed (migration, query, debug, policy change)
- Whether targeting local (`127.0.0.1:54322`) or cloud Supabase
- Any relevant table names or registry models involved
- For migrations: what schema change is needed and why

## Architecture
- Registry ORM models: `sc_tools/registry.py` (`_build_models()`)
- Alembic migrations: `sc_tools/migrations/versions/`
- Alembic config: `alembic.ini` (repo root)
- Migration env: `sc_tools/migrations/env.py`
- MCP server config: `.mcp.json` (sets `SC_TOOLS_REGISTRY_URL`)
- Data migration script: `scripts/migrate_sqlite_to_pg.py`

## Connection strings
- **Local**: `postgresql://postgres:postgres@127.0.0.1:54322/postgres`
- **Cloud**: set via `SC_TOOLS_REGISTRY_URL` env var (see `.env` or `.mcp.json`)

## Common tasks

### Start/stop local Supabase
```bash
supabase start   # spins up local Postgres + services via Docker
supabase stop    # tears down containers (data persists in Docker volumes)
supabase status  # show running services and ports
```

### Create a new Alembic migration
```bash
SC_TOOLS_REGISTRY_URL="postgresql://..." alembic revision -m "description"
# Edit the generated file in sc_tools/migrations/versions/
# Ensure revision/down_revision variables are set
SC_TOOLS_REGISTRY_URL="postgresql://..." alembic upgrade head
```

### Run migrations
```bash
SC_TOOLS_REGISTRY_URL="postgresql://postgres:postgres@127.0.0.1:54322/postgres" alembic upgrade head
```

### Migrate data from SQLite
```bash
python scripts/migrate_sqlite_to_pg.py --target "postgresql://..."
python scripts/migrate_sqlite_to_pg.py --dry-run  # preview only
```

### Direct SQL access
```bash
supabase db dump         # dump schema
psql "postgresql://postgres:postgres@127.0.0.1:54322/postgres"  # interactive shell
```

### Push schema to cloud
```bash
supabase db push         # apply local migrations to linked cloud project
supabase link --project-ref <REF>  # link to cloud project first
```

## Standards
- All schema changes go through Alembic migrations (never raw DDL in production)
- Migration files MUST have `revision` and `down_revision` module-level variables
- `from __future__ import annotations` must precede revision variables
- Migrations must be dialect-agnostic (no SQLite or PostgreSQL-specific SQL)
- Use `op.batch_alter_table` for ALTER TABLE operations (Alembic best practice)
- Test migrations on local Supabase before pushing to cloud
- Always reset sequences after bulk data loads (`setval()`)

## Verification
1. `supabase status` confirms local services are running
2. `alembic upgrade head` succeeds without errors
3. `Registry().status()` returns expected project/dataset counts
4. MCP tools (`registry_status`, `list_datasets`) work after restart
