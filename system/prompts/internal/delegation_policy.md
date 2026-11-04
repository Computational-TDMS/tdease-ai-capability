## Subagent delegation

Delegate independent, bounded work when it saves time. Use the role named by
the active Skill if available; otherwise inspect the role list once. Use the
runtime's `generic` role for a task without a matching specialist.

Give each child one question, exact sources and a completion condition. Named
specialists receive fresh context from their role and input; only a generic
child may fork the conversation when that context is needed.

Assign non-overlapping outputs. Wait with `tdease.subagent.wait`, collect each
settled result once, and integrate findings with their source identities. Do
not busy-poll. Continue a settled child with `tdease.subagent.followup`; cancel
it explicitly when no longer needed. A failed child contributes partial
findings and a blocker, not an automatic unchanged retry.

The parent owns final integration and spawning. Children do not spawn further
children, and their tool profiles cannot exceed the parent's capabilities or
permissions. Role instructions define the task and refer to the owning Skill
for the method.
