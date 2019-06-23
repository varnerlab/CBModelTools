## Expanded Varnerlab Flat File (VFF) format
The expanded Varnerlab Flat File (VFF) format is a simple text format used to work with constraint based models. It is composed of several sections, each containing different types of information:

    // Contents: Chemical reaction string
    // Records: Reaction_name,{ecnumbers::[]},reactants,products,reverse,forward
    // Records: MDH,ec:1.1.1.37::[],mal_L[c]+nad[c],h[c]+nadh[c]+oaa[c],-inf,inf
    // Records: SUCD1m,ec:3.2.1.48::ec:3.2.1.20::ec:3.2.1.3::ec:3.2.1.10::[],q10[m]+succ[m],fum[m]+q10h2[m],0,inf
    #REACTION::START --------------------------------------------- //
    ... reaction records go here ...
    #REACTION::END ----------------------------------------------- //

    // Contents: Mapping between metabolite symbol and actual name
    // Contents: and KEGG compound ID
    // Records: Metabolite_symbol=Metabolite_name::KEGG_ID
    // Records: q10[m]=Ubiquinone-10::cpd:C11378
    #METABOLITES::START ------------------------------------------ //
    ... metabolite records go here ...
    #METABOLITES::END -------------------------------------------- //

    // Contents: Rule records exported from the MAT binary file
    // Records: Reaction_name={rule} or [] if none. The rules are
    // Records: boolean statements of gene states
    // Records: SUCD1m=(x(1481) | x(1918))
    // Records: ASPte=[]
    #RULES::START ------------------------------------------------ //
    ... rules records go here ...
    #RULES::END -------------------------------------------------- //

    // Contents: Mapping between reactions in the model and genes
    // Records: Reaction_name={:: delimited gene set}::[]
    // Records: MTHFC=hsa:4522.1::hsa:80068.1::hsa:285216.1::[]
    #REACTION-GENE-MAP::START ------------------------------------ //
    ... reaction-gene-map records go here ...
    #REACTION-GENE-MAP::END -------------------------------------- //

    // Contents: Gene symbol. Order matters - add new genes to bottom
    // Records: kegg_organism_code:gene_symbol.splice
    // Records: hsa:10195.1
    #GENE-SYMBOL-ORDER::START ------------------------------------ //
    ... gene-symbol-order records go here ...
    #GENE-SYMBOL-ORDER::END -------------------------------------- //


### REACTION section
The reaction section holds chemical reaction records, which are the core of the model.
Chemical reaction records contain six fields that are comma delimited; Reaction_name,{ecnumbers::[]},reactants,products,reverse,forward.

| Field | Description |
| ----- | ----------- |
| Reaction_name | Reaction name string used in the model to denote this reaction. Must be unique
| EC numbers | [EC numbers](https://en.wikipedia.org/wiki/Enzyme_Commission_number) associated with this reaction. Estimated from KEGG using the gene linkage information. Multiple possible assignments separated with the ``::`` delimiter. If no information then ``[]``.
| Reactants | Chemical reactant string e.g., 2*A+B. Non-unity stoichiometric coefficients are affixed to metabolite by *. No spaces between entries and the ``+``.
| Products | Chemical product string e.g., 4*C+3*D. Same rules as reactant strings
| Reverse | Reaction lower bound. If reversible ``-inf`` if non-reversible 0
| Forward | Reaction upper bound. Can be either ``inf`` or 0

### METABOLITES section
The metabolites section holds records that map the metabolite symbols used in the ``REACTION`` chemical reaction strings to their actual chemical name, and an associated [KEGG metabolite ID](https://www.genome.jp/kegg/compound/).
Records in the ``METABOLITES`` section have three fields separated by ``=`` and ``::`` delimiters.

| Field | Description |
| ----- | ----------- |
| Metabolite symbol | Short-hand symbol used in the model for a metabolite name. Typically no spaces. Model compartment is indicated by a ``[*]`` suffix where `*` denotes a 1-letter compartment code e.g., ``c`` indicates cytosol. These codes are often model specific.
| Metabolite name | IUPAC name for the metabolite (most of the time).
| KEGG Metabolite ID | [KEGG compound code](https://www.genome.jp/kegg/compound/) to uniquely identify the chemical compound. If [KEGG metabolite ID](https://www.genome.jp/kegg/compound/) not provided, we estimate these using the [KEGG.jl](https://github.com/varnerlab/Kegg) package.

### RULES section
### REACTION-GENE-MAP section
### GENE-SYMBOL-ORDER section
