"""
Optional vulture whitelist for false positives (dynamic imports, CLI entry points, etc.).

Pass to vulture as a **second path** (not a flag):

    vulture . vulture_whitelist.py ...

See https://github.com/jendrikseipp/vulture#usage
"""

whitelist = [
    # Example: "SomeClass.unused_but_public_api",
]
