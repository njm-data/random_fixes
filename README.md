# random_fixes
Random pipeline fixes/band aids that I find useful.

formal_rename.py - intent: from ToxCodAn, taking the toxin_class rsem names and tallying them up (ie Toxin_extenderContigX||XXXXX_SVSP to SVSP-1) and then rename in the fasta as well. Could apply to any excel file with a long gene ID (including nontoxin NCBI names) and transitioning to a short ID based on toxin class or NCBI names.

aa_mutation_check.py - intent: compares translated and aligned sequences for differences in amino acid families.
