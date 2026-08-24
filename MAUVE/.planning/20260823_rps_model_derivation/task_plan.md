# Task plan: full derivation from shared RPS-model chat

Date: 2026-08-23

## Goal

Understand the complete linked ChatGPT conversation, identify the precise model and requested result, reconstruct the derivation from first principles with explicit assumptions and limiting cases, validate the algebra/numerics independently, and deliver a dated Markdown report and a visually checked PDF in `/Users/Igniz/Desktop/ICRAR/MAUVE`.

## Deliverables

- `20260823 <descriptive title>.md`
- `20260823 <descriptive title>.pdf`
- Supporting figure/source/build assets only when required for reproducibility

## Phases

1. **Conversation recovery and scope lock** - complete
   - Recover all user and assistant messages from the linked chat.
   - Record equations, definitions, references, assumptions, and unresolved questions verbatim or faithfully paraphrased.
   - Define what “full detailed derivation” must cover.
2. **Source and model verification** - complete
   - Check cited/source equations against primary literature or authoritative sources.
   - Separate source-backed definitions from new derivation choices.
3. **Independent derivation and computational checks** - complete
   - Derive every transformation step-by-step.
   - Check dimensions, signs, normalization, boundary/initial conditions, limiting cases, and numerical behavior.
   - Use Mathematica if available and useful; otherwise use symbolic/numerical Python with disclosed tooling.
4. **Figures and report drafting** - complete
   - Create only figures that materially clarify geometry, causal structure, or solution behavior.
   - Write the dated Markdown report with notation table, derivation, checks, interpretation, limitations, and references.
5. **PDF build and acceptance QA** - complete
   - Generate PDF, extract text, inspect metadata/page count, render every page, and visually inspect all pages.
   - Verify equations, figures, tables, links, and references are not clipped or malformed.
6. **Handoff and job log** - complete
   - Append the required concise Codex job log.
   - Report exact files, checks, and remaining limitations.

## Scope guardrails

- The linked conversation is the primary source for the requested mathematical target.
- The 2026-08-13 MAUVE RPS review is contextual evidence only; it must not replace the linked conversation.
- Do not silently invent missing equations, parameter choices, or boundary conditions.
- Clearly label exact identities, model assumptions, approximations, and optional extensions.
- Treat galaxy-history inference as probabilistic and model-dependent.

## Errors encountered

| Error | Attempt | Resolution |
|---|---:|---|
| Conversation page rendered only the ChatGPT shell, with disabled composer and no messages | 1 | Retrying after checking fresh visible state; authentication appears present. |
| Browser page-state query timed out after waiting | 2 | Next attempt will use a cheaper state check and inspect screenshot/URL before any alternative route. |
| Chrome does not support browser-native page export | 1 | Use message-specific DOM extraction in bounded chunks instead. |
| First predecessor-chat turn query timed out during load | 1 | Check title/URL/screenshot first, then retry only after visible content is present. |
| Second predecessor chat initially blank; delayed turn query timed out | 1 | Use a cheap title/screenshot check after a longer natural interval, avoiding an immediate selector query. |
| Mathematica kernel and notebook front end were unavailable | 2 | Disclosed the limitation and used a standard-library numerical verification script instead. |
| Initial PDF route preserved display equations but lost inline TeX glyphs | 2 | Rebuilt through Pandoc native MathML with the single-backslash math extension and rendered in a fresh headless-Chrome profile. |
| First numbered-page render lacked its destination directory | 1 | Created the exact QA directory and reran the non-destructive render. |

## Acceptance criteria

- Complete linked chat recovered or any inaccessible portions explicitly listed.
- Every symbol defined before use and every nontrivial algebraic step shown.
- Dimensional analysis and at least one independent symbolic/numerical check pass.
- Source claims trace to primary/authoritative material with direct links or precise document locators.
- Markdown and PDF share the same substantive content.
- Every PDF page is rendered and visually inspected with no clipping, overlap, missing glyphs, or broken equations.
