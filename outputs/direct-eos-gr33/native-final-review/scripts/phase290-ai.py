"""Phase 290: add the model versions confirmed by the author (Claude Opus 5.5, GPT-6-Astra) to the Use of AI tools section."""
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
p = root/'paper/manuscript.md'; t = p.read_bytes().decode('utf-8')
old = ("Coding agents built on OpenAI Codex models and on Anthropic Claude models, run through the Codex and Claude Code command-line tools, "
       "wrote and ran most of the analysis, stellar-structure and verification code, carried out the computations reported here, "
       "and drafted and revised manuscript text, including scientific claims and explanations. "
       "The final text of Section 4.6 and its revisions were prepared with Claude Opus 5.5 in Claude Code. "
       "GPT-6-Astra (through the Codex command-line tool), Claude Opus 5.5 and Claude Fable 5.1 independently reviewed drafts of Section 4.6 in several rounds;")
new = ("Coding agents built on Claude Opus 5.5 (Anthropic, run through Claude Code) and GPT-6-Astra (OpenAI, run through the Codex command-line tool) "
       "wrote and ran most of the analysis, stellar-structure and verification code, carried out the computations reported here, "
       "and drafted and revised manuscript text, including scientific claims and explanations. "
       "GPT-6-Astra, Claude Opus 5.5 and Claude Fable 5.1 independently reviewed drafts of Section 4.6 in several rounds;")
assert t.count(old) == 1
p.write_bytes(t.replace(old, new).encode('utf-8')); print('ok')
