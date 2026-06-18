---
name: feedback-no-modify-scripts
description: Do not modify existing scripts; create new ones instead
metadata:
  type: feedback
---

Do not modify existing scripts in this project. When new functionality is needed, create a new script (even if it means copying code from existing ones).

**Why:** User explicitly requested this to preserve the integrity of their existing pipeline scripts.

**How to apply:** Any time a task could be done by editing an existing `.py` script, instead create a new script and import from or copy the relevant functions into it.
