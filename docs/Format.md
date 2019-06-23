### Expanded Varnerlab Flat File (VFF) format
The expanded Varnerlab Flat File (VFF) format is a simple text format that can used to work with constraint based models. It is composed of several sections that each contain different types of information:

    // Contents: Mapping between reactions in the model and genes
    // Records: Reaction_name={:: delimited gene set}::[]
    // Records: MTHFC=hsa:4522.1::hsa:80068.1::hsa:285216.1::[]
    #REACTION-GENE-MAP::START ------------------------------------ //
    ... reaction-gene-map records go here ...
    #REACTION-GENE-MAP::END -------------------------------------- //

    // Contents: Rule records exported from the MAT binary file
    // Records: Reaction_name={rule} or [] if none. The rules are
    // Records: boolean statements of gene states
    // Records: SUCD1m=(x(1481) | x(1918))
    // Records: ASPte=[]
    #RULES::START ------------------------------------------------ //
    ... rules records go here ...
    #RULES::END -------------------------------------------------- //

    // Contents: Mapping between metabolite symbol and actual name
    // Contents: and KEGG compound ID
    // Records: Metabolite_symbol=Metabolite_name::KEGG_ID
    // Records: q10[m]=Ubiquinone-10::cpd:C11378
    #METABOLITES::START ------------------------------------------ //
    ... metabolite records go here ...
    #METABOLITES::END -------------------------------------------- //

    // Contents: Chemical reaction string
    // Records: Reaction_name,{ecnumbers::[]},reactants,products,reverse,forward
    // Records: MDH,ec:1.1.1.37::[],mal_L[c]+nad[c],h[c]+nadh[c]+oaa[c],-inf,inf
    // Records: SUCD1m,ec:3.2.1.48::ec:3.2.1.20::ec:3.2.1.3::ec:3.2.1.10::[],q10[m]+succ[m],fum[m]+q10h2[m],0,inf
    #REACTION::START --------------------------------------------- //
    ... reaction records go here ...
    #REACTION::END ----------------------------------------------- //

    // Contents: Gene symbol. Order matters - add new genes to bottom
    // Records: kegg_organism_code:gene_symbol.splice
    // Records: hsa:10195.1
    #GENE-SYMBOL-ORDER::START ------------------------------------ //
    ... gene-symbol-order records go here ...
    #GENE-SYMBOL-ORDER::END -------------------------------------- //
