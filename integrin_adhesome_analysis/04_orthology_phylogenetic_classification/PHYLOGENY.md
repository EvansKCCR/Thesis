# Orthology and phylogeny preparation

Follow [the step-by-step guide](../../README/ORTHOLOGY_PHYLOGENY_GUIDE.md) from the
workflow root. Phylogeny is disabled by default. Add `--phylogeny` to resume/continue;
use both `--phylogeny --phylogeny-grouping family` for unified family trees.
The audit completes before optional tree preparation/inference and remains
available if tree tools or QC fail. Use `recover --run PATH --phylogeny` to opt in
to or retry trees while reusing verified orthology, benchmark and audit outputs.
