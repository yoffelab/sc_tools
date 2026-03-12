Things to do
	sc_tools
		1. Redesign sc_tools mcp server library based on db and OOP
			1. Define abstract dataclasses such as BioData that can be extended into multiple data formats. Suggest if you have more modular objects
				1. BioImages
					1. H&E
					2. Fluoresence
					3. Multiplexed
				2. RNA-seq
					1. Bulk
					2. SingleCell
				3. SpatialSeq
					1. Bulk
					2. SingleCell
				4. Epigenomics
				5. GenomeSeq
			2. Define projects as collections (sets of these projects). Im not sure how these would go into the db. Each projects should have the capacity to work across multiple data sources.
			3. Define patient database with de-identified id (without PHI). This should hold standardized clinical metadata.
		2. Improve Deconvolution Figures
		3. Clean up documentation
			1. Move documentation related markdowns from ./ to ./docs/  (everything except for ./claude.md)
			2. Make a folder for missions (plans) to create tiered missions. Always keep and store plans. Keep copies of plan versions (V1, V1.1, V1.2, V2, etc)
				1. Tier 1 large scale overall plans (overall data flow), changes in architectural plans should go here.
				2. Tier 2 small size sub plans should go here
				3. Untiered plan should go here
		4. Add subset celltype and re-perform after cell typing
	Project Specific
		5. Robin
			1. Process additional data
			2. Point of analysis:
			   When you analyze the samples, look at the expression of Adenosine receptor genes:  ADORA2A and ADORA2B, and at the signature of the downstream signaling of the Gs-coupled receptors A2AR and A2BR: CREM, FOS, JUN, DUSP1, PDE4B

We will use this to compare with the second trial (PANTHER), where patients received a treatment that inhibits these receptors. 

This table was sent to me by Encouse (the clinician that runs PANTHER)

|   |   |   |
|---|---|---|
|Cell category|Primary receptor|What to evaluate|
|Lymphocytes (T, NK, B)|ADORA2A|Check if expression correlates with exhaustion markers; **high expression + high exhaustion = cell that was heavily suppressed**|
|Myeloid (Monocytes, DC, MDSC)|ADORA2B|Check if A2B is induced with radiation; does radiation create problem that etruma is meant to fix|
|Progenitors/stem cells|ADORA2A, ADORA2B|Adenosine can inhibit hematopoietic stem cell mobilization. Check if etruma can increase the mobilization of hematopoietic stem cells from the bone marrow.|
		1. GGO 
			1. Extended analysis:
			Thank you again for your previous analysis. I am writing to follow up and ask whether it is possible to do additional analysis on Slc2a1+ macrophages in the lung cancer patient samples.
			Our recent scRNA seq data show a significant alteration in CD8 and CD4 T cell populations in mouse PyMT breast tumors when Slc2a1+ macrophages are deficient.
			Although Slc2a1+ macrophages appear to be relatively rare in lung cancer as you mentioned before, we thought it may be worth checking whether there is any correlation between Slc2a1+ macrophages and T cells.
			So, we would like to ask if it would be possible to check the potential relationship between Slc2a1+ macrophages and CD8 T/CD4 T cells, similar to the analysis you did for Slc16a3+ neutrophils. It would also be great if this could be done in both the scRNA-seq data and the spatial tx data, if feasible.
