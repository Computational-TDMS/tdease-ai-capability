# Structure view maintenance

The shipped structure-ptm route uses offline 3Dmol. Fetching, coordinate
mapping and SASA calculation are outside this skill's scope.

Keep the actual structure, chain/residue identifiers and mapped site list in
the data contract. Label source accession/PDB, model version and coordinate
mapping. Include ligands or additional chains when they matter to the question.

Use a representation and camera that make the relevant site inspectable.
Context geometry can be muted while selected residues remain distinct.
Confidence can be shown with a legend or alongside residue values.

AlphaFold confidence is not experimental B-factor semantics. Read the source
metadata before interpreting that column. Local confidence does not resolve
inter-domain orientation; represent such uncertainty when it affects the claim.

Link 2D and 3D selections only through verified residue mapping. Preserve
camera/selection when reproducibility requires the same view. Check the WebGL
view and PNG export; SVG cannot contain the molecular canvas.
