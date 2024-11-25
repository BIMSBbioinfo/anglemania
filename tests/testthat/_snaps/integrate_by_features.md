# integrate_by_features integrates Seurat objects correctly using selected features

    Code
      SeuratObject::LayerData(se_integrated)[1:10, 1:10]
    Output
      10 x 10 sparse Matrix of class "dgCMatrix"
    Message
        [[ suppressing 10 column names 'ATGCCAGAACGACT', 'GAACCTGATGAACC', 'TGACTGGATTCTCA' ... ]]
    Output
                                                                                     
      GNG11     . .         .         .        .         .         .        .        
      CLU       . .         .         .        .         .         .        .        
      SDPR      . .         .         .        .         .         0.000000 .        
      PF4       . .         .         .        .         .         .        0.3681607
      GP9       . .         .         .        .         .         .        .        
      SPARC     . .         .         .        .         .         .        .        
      PPBP      . 1.095858 -1.474242 -1.574689 0.4642328 0.7991843 0.650469 1.4945964
      HIST1H2AC . .         .         .        1.1626601 1.2506524 1.069532 0.3354995
      CD9       . .         .         .        .         .         .        .        
      NRGN      . .         .         .        .         .         .        .        
                                   
      GNG11     .         .        
      CLU       .         .        
      SDPR      .         .        
      PF4       .         0.2570694
      GP9       .         .        
      SPARC     .         .        
      PPBP      0.7744361 0.6824838
      HIST1H2AC 1.2228468 0.8901318
      CD9       .         .        
      NRGN      .         .        

