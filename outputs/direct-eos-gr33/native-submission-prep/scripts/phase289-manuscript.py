"""Phase 289: submission edits to paper/manuscript.md (CRLF kept).

1. Drop the 'Status: Unified revised manuscript' metadata line, which the builder prints in the title block and which reads like a
   resubmission for a first submission.
2. Add a 'Use of AI tools' section before 'Data and code availability' (APS policy, June 2026: substantive AI use must be disclosed
   in the paper with tool names and versions, how the AI assisted, and how the authors directed and verified it). The author must
   confirm this text and add the model versions of earlier project stages before submitting.
"""
from pathlib import Path
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672')
p = root/'paper/manuscript.md'
t = p.read_bytes().decode('utf-8')
status = '**Status:** Unified revised manuscript\r\n'
assert t.count(status) == 1; t = t.replace(status, '')
ai = ("## Use of AI tools\r\n\r\n"
      "This work used AI agents substantively, under the author's direction and responsibility. "
      "Coding agents built on OpenAI Codex models and on Anthropic Claude models, run through the Codex and Claude Code command-line tools, "
      "wrote and ran most of the analysis, stellar-structure and verification code, carried out the computations reported here, "
      "and drafted and revised manuscript text, including scientific claims and explanations. "
      "The final text of Section 4.6 and its revisions were prepared with Claude Opus 5.5 in Claude Code. "
      "GPT-6-Astra (through the Codex command-line tool), Claude Opus 5.5 and Claude Fable 5.1 independently reviewed drafts of Section 4.6 in several rounds; "
      "their reports and the responses are archived in the repository. "
      "The author set the research questions, the acceptance criteria and the claim labels, and checked the AI output against the stored numerical records, "
      "the SHA-256-bound manifests, symbolic checks, the verification scripts listed below and these reviews.\r\n\r\n")
anchor = '## Data and code availability\r\n'
assert t.count(anchor) == 1 and 'Use of AI tools' not in t
t = t.replace(anchor, ai + anchor)
p.write_bytes(t.encode('utf-8'))
lines = t.split('\n'); print('lines', len(lines), 'bare LF', [i + 1 for i, l in enumerate(lines[:-1]) if not l.endswith('\r')])
