# Status Labels

Every scientific line in `docs/` and `lemmas/` carries a status label. The
taxonomy below fixes what each label means, so that `Proven` again marks
verified mathematics rather than being diluted across premises, open questions,
and bookkeeping.

| Label | Meaning | Examples |
| --- | --- | --- |
| `Proven` | A mathematical claim that is established -- by explicit symbolic/numeric computation in `../symbolic/`, by an independent check in `../verification/`, or by a written proof. | "The reduced scalar operator space is finite-dimensional." "`Tr(E^4) = (1/2)(Tr E^2)^2` for STF 3x3." "The seven survivors are linearly independent (rank 7)." |
| `Assumption` | A posited premise that defines the theorem domain. Neither derived nor open -- the positive theorem is *conditional* on it. Used for `A1`-`A9` (see [`assumptions-ledger.md`](assumptions-ledger.md)). | "There is no orbital-timescale internal state variable (`A4`)." |
| `Conjectural` | An open claim expected to be resolved but not yet established, or one that holds only for the exact current primitive set and is not yet a minimal-sector statement. | "This is an exact-current-set theorem candidate, not yet a physically justified minimal-sector theorem." |
| `Counterexample candidate` | A boundary/loophole model proposed as a counterexample to a stated assumption, not yet a fully closed refutation. | The dynamic-`chi` A4 escape; the non-analytic A5 activation. |
| `Note` | Provenance, file description, historical record, or bookkeeping about the documents themselves -- not a scientific claim. The whole `archive/` tree is `Note` (superseded history). | "This file is a supporting table, not an authoritative theorem statement." "This note records the assumption-drop failures." |

## Sweep applied

- `A1`-`A9`: `Conjectural` -> `Assumption` (they are premises).
- `archive/`: all `Status: Proven` -> `Status: Note` (the directory is
  quarantined, non-authoritative history).
- Active `docs/` and `lemmas/`: lines that describe the document itself
  (`This note/file/document/table/ledger/register/map/roadmap ...`) ->
  `Status: Note`. Genuine mathematical `Proven` claims are unchanged.

Lines that describe *what a lemma asks/does* (`This lemma ...`) and scope
statements (`The repository does not currently assume ...`) were left as-is;
they are borderline and can be individually reviewed, but they are not the
provenance/history lines this sweep targeted.
