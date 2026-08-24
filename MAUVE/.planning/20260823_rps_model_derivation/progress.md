# Progress: RPS-model chat derivation

## 2026-08-23

- Read Desktop and ICRAR project instructions; no MAUVE-local `AGENTS.md` exists.
- Loaded the academic-research, browser-control, PDF, and file-based planning workflows.
- Ran session catch-up for `/Users/Igniz/Desktop/ICRAR/MAUVE`; no unsynced planning context was reported.
- Inspected relevant prior MAUVE RPS review inventory and memory registry for reusable evidence/QC context.
- Created this scoped plan. Phase 1 is in progress.
- Connected to the user's signed-in Chrome session and opened the supplied ChatGPT conversation URL.
- Initial page inspection showed the ChatGPT shell without messages; the next read timed out, so conversation recovery remains in progress.
- The page then loaded successfully; visually confirmed the conversation title, the opening user premise, and the start of the assistant's comparison-table response.
- Tested a generic article selector; it returned no messages, so the next extraction will target ChatGPT's actual turn containers.
- Located four live turn containers and recovered the substantive model outline and scope-narrowing exchange.
- Attempted browser-native page export; Chrome does not support that command. Switching to bounded per-turn chunk extraction.
- Measured all four turn lengths and recovered the complete 19,603-character first assistant response in bounded chunks.
- Recovered the complete 17,036-character second assistant response; the four-turn current chat is now fully understood.
- Opened the first linked predecessor conversation; initial message-container inspection timed out while it loaded.
- Visually confirmed the first predecessor chat and measured its two turns (640-character request; 37,218-character response).
- Recovered the full 37,218-character Köppen et al. (2018) review/derivation and recorded its long-pulse, short-pulse, timescale, and 3D-extension implications.
- Opened the second predecessor chat; it is still loading and has not yet exposed message containers.
- The second predecessor chat loaded and exposed six turn containers; complete bounded extraction is next.
- Recovered the first two long sections of the second predecessor: spiral/non-axisymmetric susceptibility and the exact-edge-on orbital-dynamics problem.
- Recovered the remaining second-predecessor summary/directions. Conversation-recovery phase is substantively complete; moving to source and mathematical verification.
- Freshly checked the primary Gunn--Gott, Roediger--Brüggen, Köppen, Singh, and NFW records and recorded the exact claims that constrain the derivation.
- Completed the first-principles consistency audit of the proposed scalar criterion. Identified the arbitrary-wind torque obstruction and split the deliverable into exact vector dynamics plus an explicitly limited static susceptibility proxy.
- Completed the full independent derivation, including the beta-model/NFW environment, branch-preserving radial deprojection, disk-orientation geometry, analytic galaxy potential, hydrostatic finite-thickness gas, wind-aligned columns, exact edge-on Bessel-function chord, vector dynamics, long-pulse limits, exact arbitrary-direction short-pulse escape threshold, and observation-ready outputs.
- Attempted both Mathematica compute and notebook routes. No working kernel/front-end connection was available, so no Mathematica result is claimed; all checks were reproduced with the standard-library `verify_equations.py` script.
- Generated four explanatory figures and the reproducible figure/verification/build assets under `assets/20260823_rps_model_derivation/`.
- Wrote the final 7,541-word dated Markdown report and a tagged, 30-page A4 PDF with native MathML equations and page-number footers.
- Re-ran all numerical checks: Miyamoto--Nagai derivatives and the NFW force agree with finite differences below `6.2e-10` relative error; the exact edge-on chord agrees with direct quadrature below `2.1e-10`; point-mass extrema and impulse limiting cases agree with their analytic values.
- Rendered all 30 final PDF pages and inspected them via complete contact sheets plus full-resolution checks of the cover, contents, critical equations, figures, tables, references, and final page. No clipping, overlap, missing glyphs, or leaked TeX remained.
- Copied the checksum-matched final PDF into the requested MAUVE folder and moved intermediate PDF-build/QA scratch files to `/private/tmp/MAUVE_rps_derivation_pdf_qa_20260823`.
