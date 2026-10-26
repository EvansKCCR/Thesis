# Detection recovery in the same run — version 1.0

The recorded failure was `hmmscan` aborting at its 100,000-residue comparison
limit. `Sguin.faa` was copied from the WBPS19 **mRNA_transcripts** file, rather
than the matching protein FASTA. Its nucleotide sequences cannot be treated as
amino acids or recovered as valid protein-domain results. The completed Sman,
Shae, Sjap, Sbovis and Scura inputs are protein FASTAs and remain reusable.

From the existing workflow root in the activated Linux environment:

```bash
python adhesome_workflow.py start --resume --exclude-target Sguin --input-dir inputs/targets --targets-config targets.json --out runs/handoff --threads 4
```

The same command works after removing Sguin from `targets.json`. It continues
the existing `runs/handoff/detection/batch.json`, verifies and reuses the five
completed targets, processes the remaining configured targets, and exports the
handoff into **the same** `runs/handoff` directory. No backup or archive directory
is created by this recovery mode. Sguin's old files are left untouched; its
exclusion is recorded in the batch and final handoff manifests. Excluded targets
do not enter the candidate FASTA, domain-region import, coverage tables or target
species list.

`start --resume` applies to an unfinished detection/handoff export. The separate
`resume --handoff ...` command imports returned external context after the handoff
has completed. A sealed, completed handoff is not overwritten.

Recovery checks the completed targets' configuration, input/prepared FASTA hashes,
scan library/thresholds and classification inputs/resources before reuse. It will
not silently reuse changed evidence or omit a completed target. Pending targets
can reuse unchanged preparation and a completed domain scan with verified artifact
checksums. A failed or incomplete scan is rerun only for its pending target.

New input preflight rejects overwhelmingly nucleotide/transcript FASTAs before
output creation or external analysis. A genuine protein above HMMER's limit also
gets an explicit pre-scan error rather than a late SIGABRT. No automatic translation
or nucleotide-to-protein relabelling is performed.
