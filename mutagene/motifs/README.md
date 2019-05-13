# mutations
this is the code for Summer 2018 

## Changelog

* git repo
* I added a part to the code
* Added a library with modules
* Most recent update was to calc enrichment and mut load, which serves as a stat test
* added better import statements
  
Important Files to Run:
mutations_io.py
motifs_in_mutatgene.py

### Testing:

unittest variant_test.py - test for variant nucleotide in motif;
unittest range_change_test.py - tests to see how change in range size affects constant counts; #dysfunctional atm 
unittest overall_test.py - test for enrichment calculation;
unittest ambiguous_alt.py - test for ambiguous nuc in alt position;
unittest border_ref_test.py - test for not counting ref if at edge of seq;
unittest overlapping_mutations_test.py - test for multiple motifs in sample that share context;
unittest n_motif_test.py - test motif that contains "N";
unittest symmetrical_motif_test.py - test for overlapping sequence and rev. comp. seq + ref count on borders
unittest rev_and_forward.test.py - test for same counts with forward and reverse complimentary motifs

#### TO-DO:
- update get_enrichment to call logo creation directly
