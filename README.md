# claude_agents

A collection of agent config files for bioinformatics

## usage

### CLAUDE.md

Put this at `~/.claude/CLAUDE.md` to affect all agents in your environment.

### agents/

Set these up either globally at `~/.claude/agents` or in a project at `./.claude/agents`.

In a claude code session

```raw
/agent <worker-name> <cmd>
```

Eg

```raw
/agent bg-worker conduct PCA, including visualizing metadata factors and a scree plot, on input_data
```
