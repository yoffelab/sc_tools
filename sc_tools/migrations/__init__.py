"""Alembic migration scripts for the sc_tools registry database.

The registry uses SQLAlchemy with ``Base.metadata.create_all()`` for
initial table creation.  Alembic is used for subsequent schema migrations.

To initialise Alembic (one-time, in repo root)::

    pip install alembic
    alembic init sc_tools/migrations
    # Edit alembic.ini: script_location = sc_tools/migrations
    # Edit env.py: import sc_tools.registry models and point to your DB

To create a new migration after changing models::

    alembic revision --autogenerate -m "describe the change"
    alembic upgrade head
"""
